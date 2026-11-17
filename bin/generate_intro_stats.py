#!/usr/bin/env python3

"""
generate_intro_stats.py - Build the JSON file that feeds the statistics and "What's new"
sections of the intro page (www/intro.html).

Everything is pre-computed so the intro page never queries the database directly. Run it
from cron (daily is plenty), e.g.:

    15 3 * * * cd /path/to/gEAR/bin && ./generate_intro_stats.py

Output: www/stats_archive/intro_stats.json (written atomically)
"""

import argparse
import datetime
import html
import json
import os
import re
import subprocess
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)

import gear.db  # noqa: E402

WWW_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'www'))
DEFAULT_OUTPUT = os.path.join(WWW_DIR, 'stats_archive', 'intro_stats.json')
LANDING_DIR = os.path.join(WWW_DIR, 'landing')
PREFS_PATH = os.path.join(WWW_DIR, 'site_domain_prefs.json')

# Datasets flagged for removal shouldn't count anywhere
ACTIVE = "(d.marked_for_removal IS NULL OR d.marked_for_removal = 0)"

# Collapse the many dtype values into a few labels people recognize
DTYPE_GROUPS = {
    'single-cell-rnaseq': 'Single-cell RNA-seq',
    'single-nucleus-rnaseq': 'Single-nucleus RNA-seq',
    'spatial': 'Spatial',
    'bulk-rnaseq': 'Bulk RNA-seq',
    'microarray': 'Microarray',
    'epiviz': 'Epigenomic',
    'gosling': 'Epigenomic',
}

# Collections need at least this many datasets to show up in "What's new"
MIN_COLLECTION_SIZE = 3


def main():
    parser = argparse.ArgumentParser(description='Generate intro page statistics JSON')
    parser.add_argument('-o', '--output', default=DEFAULT_OUTPUT, help='Output JSON path')
    parser.add_argument('-n', '--recent_count', type=int, default=6, help='Items per "What\'s new" list')
    parser.add_argument('-d', '--domain_label', default=None,
                        help='Portal whose landing pages to list (default: domain_label in site_domain_prefs.json)')
    args = parser.parse_args()

    cnx = gear.db.MySQLDB().connect()
    cursor = cnx.cursor()

    result = {
        'generated': datetime.datetime.now().isoformat(timespec='seconds'),
        'totals': get_totals(cursor),
        'by_organism': get_by_organism(cursor),
        'by_dtype': get_by_dtype(cursor),
        'by_year': get_by_year(cursor),
        'recent_datasets': get_recent_datasets(cursor, args.recent_count),
        'recent_collections': get_recent_collections(cursor, args.recent_count),
        'landing_pages': get_landing_pages(args.domain_label or get_domain_label()),
    }

    cursor.close()
    cnx.close()

    tmp_path = args.output + '.tmp'
    with open(tmp_path, 'w') as fh:
        json.dump(result, fh, indent=2, default=str)
    os.replace(tmp_path, args.output)

    print(f"Wrote {args.output}")


def scalar(cursor, qry, params=None):
    cursor.execute(qry, params or ())
    row = cursor.fetchone()
    return int(row[0] or 0) if row else 0


def get_totals(cursor):
    return {
        'public_datasets': scalar(cursor, f"SELECT COUNT(*) FROM dataset d WHERE d.is_public = 1 AND {ACTIVE}"),
        'all_datasets': scalar(cursor, f"SELECT COUNT(*) FROM dataset d WHERE {ACTIVE}"),
        'organisms': scalar(cursor, f"SELECT COUNT(DISTINCT d.organism_id) FROM dataset d WHERE d.is_public = 1 AND {ACTIVE}"),
        'displays': scalar(cursor, "SELECT COUNT(*) FROM dataset_display"),
        'public_collections': scalar(cursor, "SELECT COUNT(*) FROM layout WHERE is_public = 1 OR is_domain = 1"),
        'gene_lists': scalar(cursor, "SELECT COUNT(*) FROM gene_cart"),
        'users': scalar(cursor, "SELECT COUNT(*) FROM guser"),
        'added_last_year': scalar(
            cursor,
            f"SELECT COUNT(*) FROM dataset d WHERE d.is_public = 1 AND {ACTIVE} "
            "AND d.date_added >= DATE_SUB(NOW(), INTERVAL 1 YEAR)"),
    }


def get_by_organism(cursor):
    cursor.execute(f"""
        SELECT o.label, COUNT(*)
          FROM dataset d JOIN organism o ON d.organism_id = o.id
         WHERE d.is_public = 1 AND {ACTIVE}
         GROUP BY o.label ORDER BY 2 DESC
    """)
    return [{'label': label, 'count': int(count)} for label, count in cursor.fetchall()]


def get_by_dtype(cursor):
    cursor.execute(f"SELECT d.dtype, COUNT(*) FROM dataset d WHERE d.is_public = 1 AND {ACTIVE} GROUP BY d.dtype")
    groups = {}
    for dtype, count in cursor.fetchall():
        label = DTYPE_GROUPS.get(dtype, 'Curated plots')
        groups[label] = groups.get(label, 0) + int(count)
    return [{'label': k, 'count': v} for k, v in sorted(groups.items(), key=lambda kv: -kv[1])]


def get_by_year(cursor):
    cursor.execute(f"""
        SELECT YEAR(d.date_added), COUNT(*) FROM dataset d
         WHERE d.is_public = 1 AND {ACTIVE} AND d.date_added IS NOT NULL
         GROUP BY 1 ORDER BY 1
    """)
    rows = []
    running = 0
    for year, count in cursor.fetchall():
        running += int(count)
        rows.append({'year': int(year), 'count': int(count), 'cumulative': running})
    return rows


def get_recent_datasets(cursor, limit):
    cursor.execute(f"""
        SELECT d.title, d.share_id, d.dtype, d.date_added, o.label, d.ldesc
          FROM dataset d JOIN organism o ON d.organism_id = o.id
         WHERE d.is_public = 1 AND {ACTIVE} AND d.load_status = 'completed' AND d.share_id IS NOT NULL
         ORDER BY d.date_added DESC
         LIMIT %s
    """, (limit,))
    return [{
        'title': title,
        'url': f"/p?s={share_id}",
        'dtype': DTYPE_GROUPS.get(dtype, 'Curated plots'),
        'date': date_added.date().isoformat() if date_added else None,
        'organism': organism,
        'summary': truncate(ldesc, 160),
    } for title, share_id, dtype, date_added, organism, ldesc in cursor.fetchall()]


def get_recent_collections(cursor, limit):
    # The layout table has no date, so the auto-increment id stands in for recency.
    # Requiring a few datasets filters out scratch/test collections.
    cursor.execute(f"""
        SELECT l.label, l.share_id, COUNT(DISTINCT lm.dataset_id) AS n,
               GROUP_CONCAT(DISTINCT o.label ORDER BY o.label SEPARATOR ', ')
          FROM layout l
          JOIN layout_members lm ON lm.layout_id = l.id
          JOIN dataset d ON d.id = lm.dataset_id AND d.is_public = 1 AND {ACTIVE}
          JOIN organism o ON o.id = d.organism_id
         WHERE (l.is_public = 1 OR l.is_domain = 1) AND l.share_id IS NOT NULL
         GROUP BY l.id, l.label, l.share_id
        HAVING n >= %s
         ORDER BY l.id DESC
         LIMIT %s
    """, (MIN_COLLECTION_SIZE, limit))
    return [{
        'title': label,
        'url': f"/p?l={share_id}",
        'dataset_count': int(n),
        'organisms': organisms,
    } for label, share_id, n, organisms in cursor.fetchall()]


def get_domain_label():
    """The current portal's domain_label (gear, nemo, cancergear, ...) from site_domain_prefs.json."""
    try:
        with open(PREFS_PATH) as fh:
            return json.load(fh).get('domain_label')
    except (OSError, ValueError) as e:
        print(f"WARNING: couldn't read domain_label from {PREFS_PATH}: {e}", file=sys.stderr)
        return None


def get_landing_pages(domain_label):
    """
    Each www/landing/<name>/ folder is a study landing page, described by its landing.json:

        {"portals": ["gear"], "title": "...", "summary": "...", "canonical_url": "https://..."}

    Every portal shares this directory, so only pages whose "portals" list includes the current
    domain_label (or "*") are listed. Pages without a landing.json are skipped. Recency comes
    from the folder's last git commit.
    """
    pages = []
    if not os.path.isdir(LANDING_DIR):
        return pages

    for name in os.listdir(LANDING_DIR):
        manifest_path = os.path.join(LANDING_DIR, name, 'landing.json')
        if not os.path.isfile(manifest_path):
            continue
        try:
            with open(manifest_path, encoding='utf-8') as fh:
                manifest = json.load(fh)
        except ValueError as e:
            print(f"WARNING: skipping {manifest_path}, invalid JSON: {e}", file=sys.stderr)
            continue

        portals = manifest.get('portals', [])
        if '*' not in portals and domain_label not in portals:
            continue

        thumbnail = manifest.get('thumbnail') or ''
        pages.append({'title': clean_text(manifest.get('title') or name),
                      'url': f"/landing/{name}/",
                      'summary': truncate(clean_text(manifest.get('summary', '')), 160),
                      'featured': bool(manifest.get('featured')),
                      'pmid': str(manifest.get('pmid') or ''),
                      'thumbnail': f"/landing/{name}/{thumbnail}" if thumbnail else '',
                      'date': last_changed(os.path.join(LANDING_DIR, name))})

    pages.sort(key=lambda p: p['date'], reverse=True)
    return pages


def last_changed(path):
    """
    Last git commit date for a path (a fresh checkout resets file mtimes, so they're
    useless here). Falls back to mtime outside a git checkout.
    """
    try:
        out = subprocess.run(['git', 'log', '-1', '--format=%cs', '--', path], cwd=WWW_DIR,
                             capture_output=True, text=True, timeout=10).stdout.strip()
        if out:
            return out
    except (OSError, subprocess.SubprocessError):
        pass
    return datetime.date.fromtimestamp(os.path.getmtime(path)).isoformat()


def clean_text(s):
    s = re.sub(r'<[^>]+>', ' ', s or '')
    return re.sub(r'\s+', ' ', html.unescape(s)).strip()


def truncate(s, length):
    s = clean_text(s)
    if len(s) <= length:
        return s
    return s[:length].rsplit(' ', 1)[0] + '…'


if __name__ == '__main__':
    main()
