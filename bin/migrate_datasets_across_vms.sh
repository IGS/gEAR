#!/bin/bash
#
# migrate_datasets_across_vms.sh
#
# Copies datasets (H5AD files + dataset / dataset_display rows) from a VM in one
# GCP project to a VM in another, reassigning ownership to user IDs that exist
# on the destination, and sets each dataset's newest single-gene and multigene
# displays as the owner's default displays (dataset_preference).
#
# H5AD files are copied with gcloud compute scp: down to this machine, then up
# to the destination, one dataset at a time, so only one file sits locally at once.
# A dataset's <id>.tar.gz is copied too if it exists on the source; it's optional.
#
# Run this from your own machine, where gcloud is logged in with access to BOTH
# projects.
#
# Usage:
#   ./migrate_datasets_across_vms.sh datasets.txt              # dry run (default): checks + builds SQL only
#   ./migrate_datasets_across_vms.sh datasets.txt --execute    # actually copy files and import
#
# datasets.txt: one dataset per line, whitespace-separated:
#     <dataset_id>
# Blank lines and lines starting with # are ignored. A single-dataset file works fine.
#
# Both VMs need MySQL credentials that the *mysql client* can read, not just
# mysqldump. Easiest is a [client] section in ~/.my.cnf, which both use:
#     [client]
#     user=root
#     password=whatever

set -euo pipefail

# ------------------------------------------------------------------ config ---
# Fill in the CHANGE_ME values before running (the script refuses to run otherwise).
SRC_PROJECT="CHANGE_ME-source-gcp-project"
SRC_ZONE="us-east1-b"
SRC_VM="CHANGE_ME-source-vm"
SRC_DATASET_DIR="/var/www/datasets"

DST_PROJECT="CHANGE_ME-destination-gcp-project"
DST_ZONE="us-east1-b"
DST_VM="CHANGE_ME-destination-vm"
DST_DATASET_DIR="/var/www/datasets"
DST_FILE_OWNER="www-data:www-data"     # chown applied to copied files on the destination

DB_NAME="gear_portal"
NEW_USER_ID="CHANGE_ME"	# guser.id on the destination that will own the datasets
SPACE_MARGIN_PCT=10                    # extra free space required beyond the transfer size

# plot_type values used to pick default displays for dataset_preference
SINGLE_GENE_TYPES="'scatter','tsne_static','umap_static','pca_static','tsne/umap_dynamic','bar','violin','line','svg'"
MULTIGENE_TYPES="'dotplot','heatmap','mg_violin','volcano','quadrant','mg_pca_static','mg_tsne_static','mg_umap_static'"
# -----------------------------------------------------------------------------

die()  { echo "ERROR: $*" >&2; exit 1; }
info() { echo "==> $*" >&2; }

# --- connection retries ---
# ssh and scp exit with status 255 when the connection itself fails (as opposed
# to the remote command failing), so only that status is retried. A failing
# remote command still stops the script immediately.
RETRIES=5          # total attempts per ssh/scp call
RETRY_DELAY=15     # seconds between attempts

retry() {
    local attempt=1 rc
    while true; do
        "$@" && return 0
        rc=$?
        if (( rc != 255 || attempt >= RETRIES )); then
            (( rc == 255 )) && echo "ERROR: connection failed $RETRIES times, giving up" >&2
            return $rc
        fi
        echo "WARNING: connection failed (attempt $attempt/$RETRIES); retrying in ${RETRY_DELAY}s..." >&2
        sleep "$RETRY_DELAY"
        attempt=$((attempt + 1))
    done
}

# ConnectTimeout makes a stuck connection attempt fail (and get retried) instead
# of hanging; ServerAlive* detects a connection that drops mid-command.
CONN_OPTS=(-oConnectTimeout=20 -oServerAliveInterval=30 -oServerAliveCountMax=4)
SSH_OPTS=(--quiet --ssh-flag=-T)
SCP_OPTS=(--quiet)
for o in "${CONN_OPTS[@]}"; do SSH_OPTS+=(--ssh-flag="$o"); SCP_OPTS+=(--scp-flag="$o"); done

_src_ssh() { gcloud compute ssh "$SRC_VM" --project "$SRC_PROJECT" --zone "$SRC_ZONE" "${SSH_OPTS[@]}" --command "$1"; }
_dst_ssh() { gcloud compute ssh "$DST_VM" --project "$DST_PROJECT" --zone "$DST_ZONE" "${SSH_OPTS[@]}" --command "$1"; }
_src_scp() { gcloud compute scp --project "$SRC_PROJECT" --zone "$SRC_ZONE" "${SCP_OPTS[@]}" "$@"; }
_dst_scp() { gcloud compute scp --project "$DST_PROJECT" --zone "$DST_ZONE" "${SCP_OPTS[@]}" "$@"; }
src_ssh() { retry _src_ssh "$@"; }
dst_ssh() { retry _dst_ssh "$@"; }
src_scp() { retry _src_scp "$@"; }
dst_scp() { retry _dst_scp "$@"; }

# The import feeds a file on stdin, so the redirect has to happen inside each
# attempt; otherwise a retry would start with stdin already partly consumed.
_dst_import() { _dst_ssh "mysql --default-character-set=utf8mb4 $DB_NAME" < "$1"; }

# lines in $1 that are not in $2 (both newline-separated lists)
missing_from() { comm -23 <(printf '%s\n' "$1" | sed '/^$/d' | sort -u) <(printf '%s\n' "$2" | sed '/^$/d' | sort -u); }

LIST_FILE="${1:-}"
MODE="${2:---dry-run}"
[[ -n "$LIST_FILE" && -f "$LIST_FILE" ]] || die "usage: $0 <datasets.txt> [--execute]"
[[ "$MODE" == "--dry-run" || "$MODE" == "--execute" ]] || die "unknown mode '$MODE'"

# ------------------------------------------------------------ read the list ---

for var in SRC_PROJECT SRC_ZONE SRC_VM DST_PROJECT DST_ZONE DST_VM DB_NAME NEW_USER_ID; do
    [[ "${!var}" != CHANGE_ME* ]] || die "set $var in the config section at the top of $0"
done
[[ "$NEW_USER_ID" =~ ^[0-9]+$ ]] || die "NEW_USER_ID must be numeric"

IDS=()
while read -r id _rest || [[ -n "${id:-}" ]]; do
    [[ -z "${id:-}" || "$id" == \#* ]] && continue
    [[ "$id" =~ ^[A-Za-z0-9_-]+$ ]] || die "bad dataset id: '$id'"
    IDS+=("$id")
done < "$LIST_FILE"

(( ${#IDS[@]} > 0 )) || die "no datasets found in $LIST_FILE"
dups=$(printf '%s\n' "${IDS[@]}" | sort | uniq -d)
[[ -z "$dups" ]] || die "dataset(s) listed more than once: $dups"

ALL_IDS=$(printf '%s\n' "${IDS[@]}")
IN_LIST=$(printf "'%s'," "${IDS[@]}");          IN_LIST=${IN_LIST%,}
FILES=$(printf '%s.h5ad ' "${IDS[@]}")

WORK="migration_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORK"
info "${#IDS[@]} datasets; working files in ./$WORK"

# ---------------------------------------------------------------- preflight ---
info "Checking source files..."
missing=$(src_ssh "cd '$SRC_DATASET_DIR' && for f in $FILES; do if ! test -r \"\$f\"; then echo \"\$f\"; fi; done")
[[ -z "$missing" ]] || die "H5AD files missing (or not readable by your ssh user) on source:"$'\n'"$missing"

info "Checking for optional .tar.gz files on source..."
TARBALLS=$(src_ssh "cd '$SRC_DATASET_DIR' && for f in $(printf '%s.tar.gz ' "${IDS[@]}"); do if test -r \"\$f\"; then printf '%s ' \"\$f\"; fi; done")
n_tar=$(printf '%s' "$TARBALLS" | wc -w | tr -d ' ')
info "$n_tar of ${#IDS[@]} datasets have a .tar.gz; those will be copied too"
ALL_FILES="$FILES $TARBALLS"

info "Checking source database..."
found=$(src_ssh "mysql -N -B $DB_NAME -e \"SELECT id FROM dataset WHERE id IN ($IN_LIST)\"")
notfound=$(missing_from "$ALL_IDS" "$found")
[[ -z "$notfound" ]] || die "dataset rows not found on source:"$'\n'"$notfound"

info "Checking destination for collisions..."
existing=$(dst_ssh "mysql -N -B $DB_NAME -e \"SELECT id FROM dataset WHERE id IN ($IN_LIST)\"")
[[ -z "$existing" ]] || die "datasets already exist in destination DB:"$'\n'"$existing"

existing=$(dst_ssh "for f in $ALL_FILES; do if sudo test -e '$DST_DATASET_DIR'/\"\$f\"; then echo \"\$f\"; fi; done")
[[ -z "$existing" ]] || die "files already exist in $DST_DATASET_DIR on destination:"$'\n'"$existing"

info "Checking destination user ID..."
found=$(dst_ssh "mysql -N -B $DB_NAME -e \"SELECT id FROM guser WHERE id = $NEW_USER_ID\"")
[[ -n "$found" ]] || die "user ID $NEW_USER_ID not found in destination guser table"

# Disk space. Files are staged in /tmp on the destination, then moved into the
# datasets dir. If both are on the same filesystem the move is a rename and the
# space is only needed once; otherwise each filesystem needs room for all of it.
# This machine only ever holds one file at a time, so it needs room for the largest.
info "Checking disk space..."
kb2gb() { awk -v k="$1" 'BEGIN { if (k >= 1048576) printf "%.2f GB", k / 1048576; else printf "%.1f MB", k / 1024 }'; }

sizes=$(src_ssh "cd '$SRC_DATASET_DIR' && stat -c %s $ALL_FILES")
need_kb=$(printf '%s\n' "$sizes" | awk '{ t += int(($1 + 1023) / 1024) } END { print t + 0 }')
max_kb=$(printf '%s\n' "$sizes"  | awk '{ k = int(($1 + 1023) / 1024); if (k > m) m = k } END { print m + 0 }')
need_margin_kb=$(( need_kb * (100 + SPACE_MARGIN_PCT) / 100 ))
n_files=$(printf '%s' "$ALL_FILES" | wc -w | tr -d ' ')
info "  transfer size: $(kb2gb "$need_kb") in $n_files files (largest $(kb2gb "$max_kb"))"

dst_df=$(dst_ssh "df -Pk /tmp '$DST_DATASET_DIR' | awk 'NR > 1 { print \$1, \$4 }'")
read -r tmp_dev tmp_avail <<< "$(printf '%s\n' "$dst_df" | sed -n 1p)"
read -r ds_dev  ds_avail  <<< "$(printf '%s\n' "$dst_df" | sed -n 2p)"
[[ -n "${tmp_avail:-}" && -n "${ds_avail:-}" ]] || die "could not read free space on destination"

if [[ "$tmp_dev" == "$ds_dev" ]]; then
	info "  destination free: $(kb2gb "$tmp_avail") (/tmp and datasets dir share a filesystem)"
	(( tmp_avail >= need_margin_kb )) \
        || die "not enough space on destination: need $(kb2gb "$need_margin_kb") (incl. ${SPACE_MARGIN_PCT}% margin), have $(kb2gb "$tmp_avail")"
else
	info "  destination free: /tmp $(kb2gb "$tmp_avail"), datasets dir $(kb2gb "$ds_avail") (separate filesystems)"
	(( tmp_avail >= need_margin_kb )) \
        || die "not enough space in /tmp on destination: need $(kb2gb "$need_margin_kb") (incl. ${SPACE_MARGIN_PCT}% margin), have $(kb2gb "$tmp_avail")"
	(( ds_avail >= need_margin_kb )) \
        || die "not enough space in $DST_DATASET_DIR on destination: need $(kb2gb "$need_margin_kb") (incl. ${SPACE_MARGIN_PCT}% margin), have $(kb2gb "$ds_avail")"
fi

local_avail=$(df -Pk "$WORK" | awk 'NR == 2 { print $4 }')
info "  local free: $(kb2gb "$local_avail")"
(( local_avail >= max_kb )) \
    || die "not enough local space for the largest file: need $(kb2gb "$max_kb"), have $(kb2gb "$local_avail")"

# dataset_display.id is auto-increment, so source IDs would clash with rows that
# already exist on the destination. We insert by column list without 'id' instead.
DISPLAY_COLS=$(dst_ssh "mysql -N -B -e \"SELECT GROUP_CONCAT(column_name ORDER BY ordinal_position) FROM information_schema.columns WHERE table_schema='$DB_NAME' AND table_name='dataset_display' AND column_name<>'id'\"")
[[ -n "$DISPLAY_COLS" ]] || die "could not read dataset_display columns on destination"

# ------------------------------------------------------------ dump metadata ---
info "Dumping metadata from source..."
DUMP_OPTS="--compact --no-create-info --complete-insert --default-character-set=utf8mb4 -h localhost"
src_ssh "mysqldump $DUMP_OPTS $DB_NAME dataset --where=\"id IN ($IN_LIST)\"" > "$WORK/dataset.sql"
src_ssh "mysqldump $DUMP_OPTS $DB_NAME dataset_display --where=\"dataset_id IN ($IN_LIST)\"" > "$WORK/dataset_display.sql"

# ------------------------------------------------------------ build import ---
# Rows are loaded into temporary staging tables, ownership is rewritten there,
# then copied into the real tables with foreign key checks still ON. Everything
# runs in one transaction: if any statement fails, nothing is committed.
#
# Displays are copied in source-id order, so their new IDs keep the same order.
# Because these datasets are new on the destination (checked above), every
# display row for them is one we just inserted, and MAX(id) per plot-type group
# is the newest display from the source.
#
# The final SELECTs are just a summary of what was imported and how many defaults
# were set for each plot type. Essentially a sanity-check.
{
    echo "START TRANSACTION;"
    # dataset has a FULLTEXT index, which InnoDB temporary tables can't have, so
    # copy just the columns (no indexes) instead of using LIKE.
    echo "CREATE TEMPORARY TABLE mig_dataset AS SELECT * FROM dataset WHERE 1 = 0;"
    echo "CREATE TEMPORARY TABLE mig_dataset_display LIKE dataset_display;"
    sed 's/^INSERT INTO `dataset` /INSERT INTO `mig_dataset` /' "$WORK/dataset.sql"
    sed 's/^INSERT INTO `dataset_display` /INSERT INTO `mig_dataset_display` /' "$WORK/dataset_display.sql"
    echo "UPDATE mig_dataset SET owner_id=$NEW_USER_ID;"
    echo "UPDATE mig_dataset_display SET user_id=$NEW_USER_ID;"
    echo "INSERT INTO dataset SELECT * FROM mig_dataset;"
    echo "INSERT INTO dataset_display ($DISPLAY_COLS) SELECT $DISPLAY_COLS FROM mig_dataset_display ORDER BY id;"
    cat <<SQL
INSERT INTO dataset_preference (user_id, dataset_id, display_id, is_multigene)
    SELECT ds.owner_id, dd.dataset_id, MAX(dd.id), 0
      FROM dataset_display dd JOIN dataset ds ON ds.id = dd.dataset_id
     WHERE dd.dataset_id IN ($IN_LIST) AND dd.plot_type IN ($SINGLE_GENE_TYPES)
     GROUP BY ds.owner_id, dd.dataset_id;
INSERT INTO dataset_preference (user_id, dataset_id, display_id, is_multigene)
    SELECT ds.owner_id, dd.dataset_id, MAX(dd.id), 1
      FROM dataset_display dd JOIN dataset ds ON ds.id = dd.dataset_id
     WHERE dd.dataset_id IN ($IN_LIST) AND dd.plot_type IN ($MULTIGENE_TYPES)
     GROUP BY ds.owner_id, dd.dataset_id;
SELECT 'datasets imported', COUNT(*) FROM mig_dataset
UNION ALL SELECT 'displays imported', COUNT(*) FROM mig_dataset_display
UNION ALL SELECT 'single-gene defaults set', COUNT(*) FROM dataset_preference WHERE dataset_id IN ($IN_LIST) AND is_multigene = 0
UNION ALL SELECT 'multigene defaults set', COUNT(*) FROM dataset_preference WHERE dataset_id IN ($IN_LIST) AND is_multigene = 1;
SQL
    echo "COMMIT;"
} > "$WORK/import.sql"

# safety: make sure every dumped INSERT was redirected to a staging table
if grep -qE '^INSERT INTO `dataset(_display)?` ' "$WORK/import.sql"; then
    die "some INSERTs were not redirected to staging tables; inspect $WORK/import.sql"
fi

info "Preflight OK. SQL written to $WORK/import.sql"

if [[ "$MODE" == "--dry-run" ]]; then
    info "Dry run only. Review $WORK/import.sql, then re-run with --execute."
    exit 0
fi

### Execute mode: copy files and import DB rows.

# ------------------------------------------------------------- copy files ---
# Files go first: a DB row without its H5AD is a broken dataset, while an H5AD
# without a DB row is harmless. Files land in a staging dir on the destination
# and only move into place once every one has arrived and checksums match.
STAGE="/tmp/$WORK"
LOCAL_DIR="$WORK/h5ad"
mkdir -p "$LOCAL_DIR"
dst_ssh "mkdir -p $STAGE"

n=0
for id in "${IDS[@]}"; do
    n=$((n + 1))
    info "[$n/${#IDS[@]}] $id.h5ad: downloading from source..."
    src_scp "$SRC_VM:$SRC_DATASET_DIR/$id.h5ad" "$LOCAL_DIR/"
    info "[$n/${#IDS[@]}] $id.h5ad: uploading to destination..."
    dst_scp "$LOCAL_DIR/$id.h5ad" "$DST_VM:$STAGE/"
    rm -f "$LOCAL_DIR/$id.h5ad"

    if [[ " $TARBALLS " == *" $id.tar.gz "* ]]; then
        info "[$n/${#IDS[@]}] $id.tar.gz: downloading from source..."
        src_scp "$SRC_VM:$SRC_DATASET_DIR/$id.tar.gz" "$LOCAL_DIR/"
        info "[$n/${#IDS[@]}] $id.tar.gz: uploading to destination..."
        dst_scp "$LOCAL_DIR/$id.tar.gz" "$DST_VM:$STAGE/"
        rm -f "$LOCAL_DIR/$id.tar.gz"
    fi
done
rmdir "$LOCAL_DIR"

info "Verifying checksums..."
src_sums=$(src_ssh "cd '$SRC_DATASET_DIR' && md5sum $ALL_FILES")
dst_sums=$(dst_ssh "cd $STAGE && md5sum $ALL_FILES")
printf '%s\n' "$src_sums" > "$WORK/md5_source.txt"
printf '%s\n' "$dst_sums" > "$WORK/md5_dest.txt"
[[ "$src_sums" == "$dst_sums" ]] || die "checksum mismatch (see $WORK/md5_*.txt); staged files left in $STAGE on destination, nothing imported"

info "Moving files into $DST_DATASET_DIR..."
dst_ssh "sudo chown $DST_FILE_OWNER $STAGE/* && sudo mv -n $STAGE/* '$DST_DATASET_DIR'/ && rmdir $STAGE"

# --------------------------------------------------------------- import DB ---
info "Importing metadata into destination database..."
retry _dst_import "$WORK/import.sql" | tee "$WORK/import_result.txt"

info "Done. ${#IDS[@]} datasets migrated. Logs and SQL kept in ./$WORK"