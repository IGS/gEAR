#!/opt/bin/python3

"""
Simple aggregator:
* Read a Gosling hub's groups file
* For each group, open all bigWig files listed.
* For each fixed window (e.g., 1kb), compute the mean across all bigWig files in that group. Filename is <group>_group.bw
* Write out a new bigWig file for that group with the computed means. This file can be used for that group when displaying in Gosling

Note: for large genomes this is I/O heavy; prefer running as a background job due to long runtime.
"""

import argparse
from urllib.parse import urljoin
import pyBigWig # also reads bigBed files
import requests
from pathlib import Path

# place in ../www/tracks/gosling_bigwig_groups/
OUTDIR_PATH = Path("../www/tracks/gosling_bigwig_groups/")
OUTDIR_PATH.mkdir(exist_ok=True)

def open_bws(paths):
    bws = []
    for p in paths:
        bw = pyBigWig.open(str(p))  # supports remote URLs as well.
        if not bw.isBigWig():
            raise RuntimeError(f"{p} is not a valid bigWig")
        bws.append(bw)
    return bws

def aggregate_to_bigwig(group_tracks, out_path, bin_size=50):
    """
    Aggregates multiple bigWig tracks by computing the mean signal across all tracks
    for each genomic interval, and writes the aggregated result to a new bigWig file.

    Args:
        group_tracks (list of dict): List of track dictionaries, each containing at least
            a "bigDataUrl" key pointing to a bigWig file path or URL.
        out_path (str or Path): Output file path for the aggregated bigWig file.
        bin_size (int, optional): Size of genomic bins (in bases) to use for aggregation.
            Default is 50.

    Notes:
        - Only tracks with a valid "bigDataUrl" are included in the aggregation.
        - The chromosome sizes are taken from the first bigWig file in the list.
        - For each bin, the mean value across all input bigWigs is computed.
        - Bins with a mean value of zero are skipped to reduce output file size.
        - Data is written in batches to minimize memory usage.
        - Requires the `pyBigWig` library and an `open_bws` function to open bigWig files.

    Warnings:
        - If no valid bigWig files are found in `group_tracks`, the function prints a warning and returns.
        - If a track is missing the "bigDataUrl" key, it is skipped with a warning.
    """
    bw_paths = []

    for track in group_tracks:
        data_url = track.get("bigDataUrl", None)
        if not data_url:
            print(f"WARNING: No bigDataUrl found for track '{track.get('shortLabel', '')}'; skipping.")
            continue
        bw_paths.append(data_url)

    if not bw_paths:
        print("No bigWig files found for this group; skipping aggregation.")
        return

    bws = open_bws(bw_paths)
    # use chrom sizes from first file
    chroms = bws[0].chroms()
    out = pyBigWig.open(str(out_path), "w")
    out.addHeader(list(chroms.items()))

    # accumulate entries per chrom and write in batches
    for chrom, chrom_len in chroms.items():
        starts = []
        ends = []
        values = []
        for start in range(0, chrom_len, bin_size):
            end = min(start + bin_size, chrom_len)
            # get mean for each bw for the interval
            vals = []
            for bw in bws:
                # bw.stats returns list with one element (or None)
                m = bw.stats(chrom, start, end, type="mean")[0]
                vals.append(0.0 if m is None else float(m))
            mean_val = sum(vals) / len(vals)
            # optional: skip zeros to reduce size
            if mean_val != 0.0:
                starts.append(start)
                ends.append(end)
                values.append(mean_val)
            # flush periodically to keep memory small
            if len(starts) >= 100000:
                out.addEntries([chrom]*len(starts), starts, ends=ends, values=values)
                starts = []; ends = []; values = []

        if starts:
            out.addEntries([chrom]*len(starts), starts, ends=ends, values=values)

    out.close()
    for bw in bws:
        bw.close()

def fetch_trackdb_and_groups_info(genomes_txt, assembly) -> dict:
    """Extract 'trackDb' and 'groups' URLs for an assembly from genomes_txt.

    Looks for a "genome <assembly>" line, then reads the immediate following
    "trackDb" and optional "groups" lines. Returns a dict with keys
    "trackDb" and "groups" (empty string if not found).

    NOTE: These can be relative paths to the genomes.txt location; caller
    must resolve them if needed.
    """
    urls = {"trackDb": "", "groups": ""}

    for line in genomes_txt.splitlines():
        if line.startswith(f"genome {assembly}"):
            # The next line should contain the trackDb URL
            next_line_index = genomes_txt.splitlines().index(line) + 1
            if next_line_index < len(genomes_txt.splitlines()):
                next_line = genomes_txt.splitlines()[next_line_index]
                if next_line.startswith("trackDb"):
                    urls["trackDb"] = next_line.split(" ")[1]
                    # Now look for groups line (next line)
                    next_next_line_index = next_line_index + 1
                    if next_next_line_index < len(genomes_txt.splitlines()):
                        next_next_line = genomes_txt.splitlines()[next_next_line_index]
                        if next_next_line.startswith("groups"):
                            urls["groups"] = next_next_line.split(" ")[1]
                    break
    return urls

def parse_groups_from_groupsdb(groupsdb_txt) -> dict:
    """Parse groups info from a UCSC groups.txt content.

    Returns a dict with group names as keys and their descriptions as values.

    Example group:

    name ATAC
    label ATAC-seq
    defaultIsClosed 0
    """
    groups = {}
    current_group = None
    for line in groupsdb_txt.splitlines():
        if line.startswith("name"):
            current_group = line.split(" ")[1]
            groups[current_group] = {}
        elif current_group:
            if line.startswith("label"):
                groups[current_group]["label"] = line.split(" ")[1]
            elif line.startswith("defaultIsClosed"):
                groups[current_group]["defaultIsClosed"] = line.split(" ")[1]
    return groups

def parse_tracks_from_trackdb(trackdb_txt, trackdb_url) -> list:
    """Parse track names from a UCSC trackDb.txt content.

    Returns a list of dicts with track information.

    Example track:
    track P1HC_ATAC_1
    bigDataUrl P1HC_ATAC_1.bigwig
    shortLabel ATAC-seq 1st replicate
    longLabel ATAC-seq 1st replicate
    group ATAC
    color 31,119,180
    autoscale on
    visibility dense
    type bigWig
    """

    tracks = [] # List of dicts
    current_track = {}
    for line in trackdb_txt.splitlines():
        if line.startswith("track"):
            if current_track:
                tracks.append(current_track)
            current_track = {"name": line.split(" ")[1]}
        elif current_track:
            if line.startswith("bigDataUrl"):
                current_track["bigDataUrl"] = line.split(" ")[1]
                # If not a URL, make it one by replacing the "trackDb.txt" part of trackdb_url
                if not current_track["bigDataUrl"].startswith("http://") \
                    and not current_track["bigDataUrl"].startswith("https://"):
                    current_track["bigDataUrl"] = urljoin(trackdb_url, current_track['bigDataUrl'])
            elif line.startswith("group"):
                current_track["group"] = line.split(" ")[1]
            elif line.startswith("type"):
                current_track["type"] = line.split(" ")[1]
                if current_track["type"] != "bigWig":
                    # Not a bigWig track; skip
                    current_track = {}
                    continue
    if current_track:
        tracks.append(current_track)
    return tracks


def main():

    parser = argparse.ArgumentParser(description="Aggregate BigWig files from a Gosling hub's groups file.")
    parser.add_argument("--hub_url", required=True, help="URL to the hub.txt file of the Gosling hub.")
    parser.add_argument("--assembly", required=True, help="Genome assembly (e.g., hg38, mm10).")
    parser.add_argument("--bin-size", type=int, default=50)

    args = parser.parse_args()
    hub_url = args.hub_url
    assembly = args.assembly

    # Cut off name of hub_url (hub.txt). This will be used to build more paths
    base_url = hub_url.rsplit('/', 1)[0]  # Get base URL of hub.txt

    # Look for a genomes.txt file to matching the assembly to get the "trackDb" file
    # This file is the same as the one required by the UCSC Genome Browser
    genomes_url = f"{base_url}/genomes.txt"
    print("Fetching genomes.txt...")
    try:
        genomes_response = requests.get(genomes_url)
        genomes_response.raise_for_status()
    except requests.RequestException as e:
        print(f"Error fetching genomes.txt from {genomes_url}: {e}")
        return

    print(f"Looking for trackDb and groups info for assembly {assembly}...")
    urls = fetch_trackdb_and_groups_info(genomes_response.text, assembly)
    if not urls:
        print(f"trackDb URL not found for assembly {assembly} in genomes.txt")
        return

    # If trackDb and groups URLs are relative to genomes.txt, resolve them
    trackdb_url = urls["trackDb"]
    groups_url = urls["groups"]

    if not trackdb_url.startswith("http://") and not trackdb_url.startswith("https://"):
        trackdb_url = f"{base_url}/{trackdb_url}"
    if groups_url and not groups_url.startswith("http://") and not groups_url.startswith("https://"):
        groups_url = f"{base_url}/{groups_url}"

    # Fetch tracks
    print(f"Fetching trackDb from {trackdb_url}...")
    try:
        trackdb_response = requests.get(trackdb_url)
        trackdb_response.raise_for_status()
    except requests.RequestException as e:
        print(f"Error fetching trackDb from {trackdb_url}: {e}")
        return
    tracks = parse_tracks_from_trackdb(trackdb_response.text, trackdb_url)

    # Attempt to fetch groups info (non-fatal)
    groups = {}
    if not groups_url:
        print("INFO: No groups URL provided; Exiting")
        return

    print(f"Fetching groups info from {groups_url}...")
    try:
        groups_response = requests.get(groups_url)
        groups_response.raise_for_status()
        # Parse groups info as needed; here we just store the raw text
        groups = parse_groups_from_groupsdb(groups_response.text)
    except requests.RequestException:
        print("INFO: Could not retrieve groups info; Exiting.")
        return  # Ignore errors in fetching groups

    print(f"Found {len(groups)} groups. Aggregating bigWig files for each group...")
    for group in groups:
        # TODO: Aggregate by group.  May require extra packages. Can cache (maybe upon upload even).
        # Get tracks associated with this group
        group_tracks = [track for track in tracks if track.get("group", "") == group]
        if not group_tracks:
            continue

        out_path = OUTDIR_PATH / f"{group}_group.bw"
        if out_path.is_file():
            print(f"Output file {out_path} already exists; skipping aggregation for group '{group}'.")
            continue

        print(f"Aggregating group '{group}' with {len(group_tracks)} tracks to {out_path}...")
        aggregate_to_bigwig(group_tracks, out_path, args.bin_size)


if __name__ == "__main__":
    main()