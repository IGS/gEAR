#!/opt/bin/python3

"""
Migrate saved Epiviz curations from the geardb dataset_display table to our Gosling equivalent in the same table.
We will have to do the following:
a) build the hub.txt file
b) build the genomes.txt file
c) build the tracksdb.txt file
d) build the groups.txt file
e) build the directory structure to host the above files
f) (optional) migrate the "measurement" datasourceid files to a new area.

Example epiviz display:

| id  | dataset_id                           | user_id | label  | plot_type | plotly_config |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
+-----+--------------------------------------+---------+--------+-----------+---------------+
| 885 | f0392d7e-86a8-fbf8-71b8-920e4a97877b |     204 | Epiviz | epiviz   | <see below>   |

{
    "genome": "mm10",
    "fullviewer": "http://epiviz.umgear.org/",
    "dataserver": "/epivizfile",
    "tracks": {
        "EPIVIZ-GENES-TRACK": [
            {
                "measurements": [
                    {
                        "id": "mm10",
                        "name": "mm10",
                        "type": "range",
                        "datasourceId": "https://obj.umiacs.umd.edu/genomes/mm10/mm10.txt.gz.tbi",
                        "datasourceGroup": "mm10",
                        "dataprovider": "fileapi",
                        "formula": null,
                        "defaultChartType": "GenesTrack",
                        "annotation": null,
                        "minValue": "default",
                        "maxValue": "default",
                        "metadata": [
                            "geneid",
                            "exons_start",
                            "exons_end",
                            "gene"
                        ]
                    }
                ],
                "settings": {
                    "title": "",
                    "marginTop": 25,
                    "marginBottom": 23,
                    "marginLeft": 20,
                    "marginRight": 10
                },
                "colors": [
                    "#f9a65a",
                    "#599ad3",
                    "#79c36a",
                    "#f1595f",
                    "#727272",
                    "#cd7058",
                    "#d77fb3"
                ]
            }
        ],
        "EPIVIZ-MULTISTACKED-LINE-TRACK": [
            {
                "measurements": [
                    {
                        "id": "c59ce961-9a55-4a64-a1cf-6ffc627de0b3",
                        "name": "mutant_Mes3",
                        "type": "feature",
                        "datasourceId": "/home/jorvis/git/gEAR/www/datasets_epigenetic/c59ce961-9a55-4a64-a1cf-6ffc627de0b3.bw",
                        "datasourceGroup": "c59ce961-9a55-4a64-a1cf-6ffc627de0b3",
                        "dataprovider": "fileapi",
                        "formula": null,
                        "defaultChartType": "lineTrack",
                        "annotation": null,
                        "minValue": "default",
                        "maxValue": "default",
                        "metadata": []
                    },
                    {
                        "id": "36bf69f8-de95-337c-1877-53e0cfc9ba12",
                        "name": "mutant_Mes1",
                        "type": "feature",
                        "datasourceId": "/home/jorvis/git/gEAR/www/datasets_epigenetic/36bf69f8-de95-337c-1877-53e0cfc9ba12.bw",
                        "datasourceGroup": "36bf69f8-de95-337c-1877-53e0cfc9ba12",
                        "dataprovider": "fileapi",
                        "formula": null,
                        "defaultChartType": "lineTrack",
                        "annotation": null,
                        "minValue": "default",
                        "maxValue": "default",
                        "metadata": []
                    },
                    ...
                ],
                "settings": {
                    "title": "",
                    "marginTop": 25,
                    "marginBottom": 23,
                    "marginLeft": 20,
                    "marginRight": 10,
                    "step": 1,
                    "showLines": true,
                    "showFill": true,
                    "showErrorBars": true,
                    "pointRadius": 1,
                    "lineThickness": 1,
                    "yMin": "default",
                    "yMax": "default",
                    "abLine": "default",
                    "showPoints": false
                },
                "colors": [
                    "#1f77b4",
                    "#ff7f0e",
                    "#2ca02c",
                    "#d62728",
                    "#9467bd",
                    "#8c564b",
                    "#e377c2",
                    "#7f7f7f",
                    "#bcbd22",
                    "#17becf"
                ]
            }
        ]
    }
}

Example of a gosling spec
| id   | dataset_id                           | user_id | label        | plot_type | plotly_config |
+------+--------------------------------------+---------+--------------+-----------+---------------+
| 3460 | cc6cb9ec-682c-49e0-986e-80cb80f2c582 |     662 | Gosling Test | gosling   | <see below>   |

{
    "hubUrl":"https://umgear.org/tracks/Litao_bigwigs/hub.txt",
    "assembly":"mm10",
    "goslingSpecDirUrl":"https://umgear.org/tracks/Litao_bigwigs", (probably not used)
    "gene_symbol": "Pou4f3"
}

"""

import json
from pathlib import Path
import sys
import subprocess
import shlex

from Bio import bgzf
import colorcet as cc

lib_path = Path(__file__).resolve().parent.parent / 'lib'
sys.path.append(lib_path.as_posix())
import geardb

BIGWIG_EXTENSIONS = [".bw", ".bigwig"]
BIGBED_EXTENSIONS = [".bb", ".bigbed"]

SUPPORTED_ASSEMBLIES = ["danRer10", "galGal6", "hg19", "hg38", "mm10", "rn6"]

HEX_COLORS = cc.glasbey_dark


# Set to True for local testing (localhost:8080), False for production
TEST = True
# Set to True to not actually perform any file operations or DB updates
DRY_RUN = True

def convert_hex_to_track_rgb(hex: str):
    """
    Convert a hex value into the tracksdb.txt RGB "color" values (i.e. 31,119,180)
    """

    hex = hex.lstrip('#')
    lv = len(hex)
    return ','.join(str(int(hex[i:i + lv // 3], 16)) for i in range(0, lv, lv // 3))


def bigbed_to_bed(bigbed_file, output_dir) -> bool:
    """
    Convert a bigBed file to a bed file using the UCSC tool bigBedToBed.
    Next, bgzip the file and tabix the bgzipped file.
    The bed file will be saved in the output_dir
    """

    if DRY_RUN:
        print(f"[DRY RUN] Would convert {bigbed_file} to BED in {output_dir}")
        return True

    if not any(bigbed_file.endswith(ext) for ext in BIGBED_EXTENSIONS):
        print(f"File {bigbed_file} is not a bigBed file.")
        return False

    bedbed_path = Path(bigbed_file)
    bed_path = bedbed_path.with_suffix('.bed')
    bed_path = Path(output_dir) / bed_path.name

    try:
        subprocess.run(["bigBedToBed", bigbed_file, bed_path.as_posix()], check=True)
        print(f"Converted {bigbed_file} to {bed_path.as_posix()}.")

        gz_path = bed_path.with_suffix(bed_path.suffix + '.gz')
        with bed_path.open('rb') as f_in, bgzf.BgzfWriter(gz_path) as f_out:
            # write the entirety of the input
            f_out.write(f_in.read())
        create_tabix_indexed_file(gz_path, file_type="bed")

        return True
    except subprocess.CalledProcessError as e:
        print(f"Error converting {bigbed_file} to bed: {e}")
        return False

def build_genomes_file(genomes_path, genome) -> bool:
    """ Example
    genome mm10
    trackDb mm10/trackDb.txt
    groups groups.txt
    """

    if DRY_RUN:
        print(f"[DRY RUN] Would write genomes.txt at {genomes_path} for genome {genome}")
        return True

    with open(genomes_path, "w") as f_out:
        f_out.write(f"genome {genome}\n")
        f_out.write(f"trackDb {genome}/trackDb.txt\n")
        f_out.write("groups groups.txt\n")
    return True

def build_groups_file(groups_path, groups) -> bool:
    """ Example
    name ATAC
    label ATAC-seq
    defaultIsClosed 0
    """
    if DRY_RUN:
        print(f"[DRY RUN] Would write groups.txt at {groups_path} for groups: {groups}")
        return True

    with open(groups_path, "w") as f_out:
        for group in list(groups):
            f_out.write(f"name {group}\n")
            f_out.write(f"label {group}\n")
            f_out.write(f"defaultIsClosed 0\n")
            f_out.write("\n")
    return True

def build_trackdb_file(trackdb_path: Path, tracks: list, groups: set = set()) -> bool:
    """ Example
    track P1HC_H3k4m1_2
    bigDataUrl P1HC_H3k4m1_2.bigwig
    shortLabel H3k4m1 2nd replicate
    longLabel Histone protein H3 (Lysine 4) monomethylation 2nd replicate
    group H3K4m1
    color 44,160,44
    autoscale on
    visibility dense
    type bigWig
    """

    if DRY_RUN:
        print(f"[DRY RUN] Would write trackDb.txt at {trackdb_path} for tracks: {[t['id'] for t in tracks]}")
        for track in tracks:
            print(f"[DRY RUN] Would process track: {track['id']} ({track['datasourceId']})")
            if any(track['datasourceId'].endswith(ext) for ext in BIGBED_EXTENSIONS):
                print(f"[DRY RUN] Would convert bigBed to bed: {track['datasourceId']}")
            print(f"[DRY RUN] Would copy file to {trackdb_path.parent}")
        return True

    trackdb_root = trackdb_path.parent.as_posix()

    groups_list = list(groups)

    # For simplicity, we will create a single track entry for the dataset.
    with open(trackdb_path, "w") as f_out:
        for track in tracks:
            # Get track type based on file extension
            if any(track['datasourceId'].endswith(ext) for ext in BIGWIG_EXTENSIONS):
                track['type'] = "bigWig"
            elif any(track['datasourceId'].endswith(ext) for ext in BIGBED_EXTENSIONS):
                track['type'] = "bigBed"
                are_files_made = bigbed_to_bed(track['datasourceId'], trackdb_root)
                if not are_files_made:
                    print(f"Failed to convert bigBed file {track['datasourceId']} to bed format. Skipping this track.")
                    continue
            else:
                print(f"Unsupported file type for datasourceId: {track['datasourceId']}. Skipping this track")
                continue

            # Copy file over to trackdb_root
            # For BigBed, we keep the BigBed path for exporting to UCSC,
            # but in the gosling API we will search for the Bed equivalent
            data_filename = Path(track['datasourceId']).name
            new_data_path = Path(trackdb_root) / data_filename
            Path(track['datasourceId']).replace(new_data_path)
            # Update datasourceId to just the filename (relative to tracksDb.txt)
            track['datasourceId'] = data_filename

            f_out.write(f"track {track['id']}\n")
            f_out.write(f"bigDataUrl {track['datasourceId']}\n")
            f_out.write(f"shortLabel {track['name']}\n")
            f_out.write(f"longLabel {track['name']}\n")

            if "datasourceGroup" in track:
                f_out.write(f"group {track['datasourceGroup']}\n")
                # Assign a color based on the track index
                try:
                    group_index = groups_list.index(track['datasourceGroup'])
                    if group_index:
                        hex_color = HEX_COLORS[group_index]
                        track_color = convert_hex_to_track_rgb(hex_color)
                        f_out.write(f"color {track_color}\n")
                except ValueError:
                    pass

            #f_out.write(f"color {track['color']}\n")
            f_out.write("autoscale on\n")
            f_out.write("visibility dense\n")
            f_out.write(f"type {track['type']}\n")
            f_out.write("\n")
    return True

def build_hub_file(hub_path, dataset_id) -> bool:
    """ Example
    hub Litao_bigwigs
    shortLabel Litao bigwigs
    longLabel Bigwig files for sandboxing by Litao Tao
    genomesFile genomes.txt
    email sadkins@som.umaryland.edu
    """

    if DRY_RUN:
        print(f"[DRY RUN] Would write hub.txt at {hub_path} for dataset {dataset_id}")
        return True

    dataset = geardb.get_dataset_by_id(dataset_id)
    if not dataset:
        print(f"Dataset with ID {dataset_id} not found.")
        return False

    user = geardb.get_user_by_id(dataset.user_id)
    if not user:
        print(f"User with ID {dataset.user_id} not found.")
        return False

    title = dataset.title if dataset.title else dataset.id

    with open(hub_path, "w") as f_out:
        f_out.write(f"hub {dataset.id.replace(' ', '_')}\n")
        f_out.write(f"shortLabel {dataset.id}\n")
        f_out.write(f"longLabel {title}\n")
        f_out.write("genomesFile genomes.txt\n")
        f_out.write(f"email {user.email}\n")

    return True

def convert_epiviz_to_gosling(display, e_config) -> dict | None:
    """
    Convert an Epiviz configuration dictionary to a Gosling specification dictionary.
    This is a stub function and should be implemented based on the specific conversion logic.

    {
    "hubUrl":"https://umgear.org/tracks/Litao_bigwigs/hub.txt",
    "assembly":"mm10",
    "goslingSpecDirUrl":"https://umgear.org/tracks/Litao_bigwigs", (probably not used)
    "gene_symbol": "Pou4f3"
    }

    """

    genome: str = e_config["genome"]

    # Could use just display_id but I like having the dataset ID in there for easy cross-referencing
    track_location = f"{display.dataset_id}_epiviz_{display.id}"

    # Should already be made when adding the various tracks metadata files
    tracks_path = Path(__file__).resolve().parent.parent / 'www' / 'tracks' / track_location

    if TEST:
        url_root = f"http://localhost:8080/tracks/{track_location}"
    else:
        url_root = get_domain_url() + f"/tracks/{track_location}"

    hub_path = tracks_path / "hub.txt"
    if not hub_path.is_file():
        if DRY_RUN:
            print(f"[DRY RUN] Hub file {hub_path.as_posix()} does not exist but should exist. Faking Gosling spec.")
            return {
                "hubUrl": url_root + "/hub.txt",
                "assembly": genome,
                "gene_symbol": "Pou4f3", # Placeholder
            }

        print(f"Hub file {hub_path.as_posix()} does not exist but should exist. Cannot create Gosling spec.")
        return None

    hub_url = url_root + "/hub.txt"

    gosling_spec = {
        "hubUrl": hub_url,
        "assembly": genome,
        "gene_symbol": "Pou4f3", # Placeholder
        # We don't need goslineSpecUrl
    }

    return gosling_spec

def create_tabix_indexed_file(bgzipped_file, file_type="bed") -> None:
    """
    Tabix-indexes a genomic file.
    bgzipped_file: Path to the bgzipped input file.
    file_type: Type of the file (e.g., "bed", "vcf", "gff").
    """

    if DRY_RUN:
        print(f"[DRY RUN] Would tabix-index {bgzipped_file} as {file_type}")
        return

    tabix_cmd = f"tabix -s 1 -b 2 -e 3 -p {file_type} {bgzipped_file}"
    subprocess.run(shlex.split(tabix_cmd), check=True)

def create_ucsc_style_track_files(display, e_config) -> bool:
    """
    Create UCSC-style track files (hub.txt, genomes.txt, tracksdb.txt, groups.txt)
    based on the Epiviz configuration.
    This is a stub function and should be implemented based on the specific file formats.
    """
    genome: str = e_config["genome"]
    tracks = e_config["tracks"]

    # Could use just display_id but I like having the dataset ID in there for easy cross-referencing
    track_location = f"{display.dataset_id}_epiviz_{display.id}"

    tracks_path = Path(__file__).resolve().parent.parent / 'www' / 'tracks' / track_location
    genomes_dir = tracks_path / genome
    genomes_dir.mkdir(parents=True, exist_ok=True)

    hub_path = tracks_path / "hub.txt"
    genomes_path = tracks_path / "genomes.txt"
    groups_path = tracks_path / "groups.txt"
    trackdb_path = genomes_dir / "trackDb.txt"

    if DRY_RUN:
        print(f"[DRY RUN] Would create UCSC-style track files for display {display.id} at {tracks_path}")

    # build hub.txt
    is_success = build_hub_file(hub_path, display.dataset_id)
    if not is_success:
        print(f"Failed to build hub file for display ID {display.id}")
        return False

    # build genomes.txt
    is_success = build_genomes_file(genomes_path, genome)
    if not is_success:
        print(f"Failed to build genomes file for display ID {display.id}")
        return False

    # build groups.txt
    groups = set()
    for t in tracks:
        if "datasourceGroup" in t:
            groups.add(t.get("datasourceGroup"))

    is_success = build_groups_file(groups_path, groups)
    if not is_success:
        print(f"Failed to build groups file for display ID {display.id}")
        return False

    # build trackDb.txt
    is_success = build_trackdb_file(trackdb_path, tracks, groups)
    if not is_success:
        print(f"Failed to build trackDb file for display ID {display.id}")
        return False

    return True


def get_displays(cursor):
    """
    Populates the dataset displays attribute, a list of DatasetDisplay objects
    related to this dataset.
    """

    displays = []

    qry = """
            SELECT id, user_id, dataset_id, label, plot_type, plotly_config
            FROM dataset_display
            WHERE plot_type = 'epiviz'
    """
    cursor.execute(qry)
    displays.clear()

    rows = cursor.fetchall()
    if rows is None or len(rows) == 0:
        cursor.close()
        return displays

    for display_id, user_id, dataset_id, label, plot_type, plotly_config in rows:
        display = geardb.DatasetDisplay(
            id=display_id,
            dataset_id=dataset_id,
            user_id=user_id,
            label=label,
            plot_type=plot_type,
            plotly_config=plotly_config,
        )

        displays.append(display)

    cursor.close()

    return displays

def get_domain_url():
    json_path = Path(__file__).resolve().parent.parent / 'www' / 'site_domain_prefs.json'
    json_dict = json.load(open(json_path, 'r'))
    return json_dict.get("domain_url")

def insert_gosling_display(cursor, dataset_id, user_id, gosling_config) -> int | None:
    """
    Insert a new dataset_display entry with the gosling configuration.
    Returns the new display ID or None on failure.
    """

    if DRY_RUN:
        print(f"[DRY RUN] Would insert Gosling display for dataset {dataset_id}, user {user_id}")
        return 999999  # Dummy ID for dry run

    qry = """
        INSERT INTO dataset_display (dataset_id, user_id, label, plot_type, plotly_config)
        VALUES (%s, %s, %s, %s, %s)
        RETURNING id
    """
    label = f"Gosling Migration of {dataset_id}"
    plot_type = "gosling"
    plotly_config_json = json.dumps(gosling_config)

    try:
        cursor.execute(qry, (dataset_id, user_id, label, plot_type, plotly_config_json))
        new_id = cursor.fetchone()[0]
        return new_id
    except Exception as e:
        print(f"Error inserting gosling display: {e}")
        return None

def parse_epiviz_plotly_config(config: str) -> dict | None:
    """
    Parse the Epiviz plotly_config JSON to extract relevant information for migration.
    """

    try:
        config_data = json.loads(config)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        return None

    # Extract genome and tracks information
    genome = config_data.get("genome")
    if not genome:
        print("No genome specified in Epiviz configuration.")
        return None

    tracks = config_data.get("tracks", {})
    if not tracks:
        print("No tracks found in Epiviz configuration.")
        return None

    epiviz_config = {"genome": genome, "tracks": []}

    # tracks is a dict
    line_tracks = tracks.get("EPIVIZ-MULTISTACKED-LINE-TRACK", [])
    if line_tracks:
        measurements = line_tracks[0].get("measurements", [])
        if measurements:
            epiviz_config["tracks"] = measurements
            return epiviz_config
        else:
            print("No measurements found in line tracks.")
    else:
        print("No line tracks found in Epiviz configuration.")

    return None

def update_dataset_dtype(cursor, dataset_id):
    if DRY_RUN:
        print(f"[DRY RUN] Would update dataset dtype to 'gosling' for dataset {dataset_id}")
        return None

    # Update the dataset dtype from "epiviz" to "gosling"
    update_dtype_query = """
        UPDATE dataset
        SET dtype = 'gosling'
        WHERE id = %s
    """
    cursor.execute(update_dtype_query, (dataset_id,))

def update_dataset_preferences(cursor, display_id, new_display_id):

    if DRY_RUN:
        print(f"[DRY RUN] Would update dataset_preference entries from display ID {display_id} to {new_display_id}")
        return None

    # Update any dataset_preference entries with this display_id to use the new display_id entry
    update_preference_query = """
        UPDATE dataset_preference
        SET display_id = %s
        WHERE display_id = %s
    """
    cursor.execute(update_preference_query, (new_display_id, display_id))

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    try:
        displays = get_displays(cursor)
        if not displays:
            print("No Epiviz displays found.")
            return

        for display in displays:

            print(f"--- Processing display ID {display.id}...")

            epiviz_config = parse_epiviz_plotly_config(display.plotly_config)
            if not epiviz_config:
                print(f"Skipping display ID {display.id} due to configuration issues.")
                continue


            genome = epiviz_config["genome"]
            if genome not in SUPPORTED_ASSEMBLIES:
                print(f"Skipping display ID {display.id} due to unsupported genome assembly: {genome}.")
                continue

            is_success = create_ucsc_style_track_files(display, epiviz_config)
            if not is_success:
                print(f"Skipping display ID {display.id} due to track file creation issues.")
                continue

            gosling_config = convert_epiviz_to_gosling(display, epiviz_config)
            if not gosling_config:
                print(f"Skipping display ID {display.id} due to Gosling conversion issues.")
                continue

            # Insert a new dataset_display with the gosling config
            new_display_id = insert_gosling_display(cursor, display.dataset_id, display.user_id, gosling_config)

            update_dataset_preferences(cursor, display.id, new_display_id)

            update_dataset_dtype(cursor, display.dataset_id)

            print(f"Successfully parsed Epiviz display ID {display.id} for genome {epiviz_config['genome']}.")
    except Exception as e:
        print(f"Error retrieving displays: {e}")
    finally:
        cnx.close()



if __name__ == "__main__":
    main()

