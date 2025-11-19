"""
Create a MySQL DB dump of only datasets referenced by layout_displays for given layout ids.

Usage:
    python create_test_mysql_dump.py --layout-ids 1 2 3 --dump-file /tmp/mini-gear.sql

Notes:
 - This script uses the project's DB config (ServerConfig) and lib.geardb.Connection to connect.
 - It is destructive. Use --dry-run first, then --yes to execute without interactive prompt.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import typing
from typing import Iterable

lib_path = Path(__file__).resolve().parent.parent / 'lib'
sys.path.append(lib_path.as_posix())

from gear.serverconfig import ServerConfig  # project config
from geardb import Connection  # use project's DB connection wrapper

if typing.TYPE_CHECKING:
    import configparser

CHUNK = 1000


def chunked_iter(items: Iterable, size: int = CHUNK):
    it = iter(items)
    while True:
        chunk = []
        try:
            for _ in range(size):
                chunk.append(next(it))
        except StopIteration:
            if chunk:
                yield chunk
            break
        yield chunk

def sql_quote(val: str) -> str:
    """
    Quotes a string value for use in SQL statements.

    This function wraps the input string in single quotes and escapes any existing single quotes
    by doubling them, following SQL string literal conventions.

    Args:
        val (str): The string value to be quoted.

    Returns:
        str: The SQL-quoted string.
    """
    return "'" + str(val).replace("'", "''") + "'"

def get_display_and_dataset_ids(conn: Connection, layout_ids: list):
    cur = conn.get_cursor()
    sql = "SELECT DISTINCT display_id FROM layout_displays WHERE layout_id IN ({})".format(
        ",".join("%s" for _ in layout_ids)
    )
    cur.execute(sql, tuple(layout_ids))
    rows = cur.fetchall() or []
    display_ids = {r[0] for r in rows}

    if not display_ids:
        cur.close()
        return set()

    sql2 = "SELECT DISTINCT dataset_id FROM dataset_display WHERE id IN ({})".format(
        ",".join("%s" for _ in display_ids)
    )
    cur.execute(sql2, tuple(display_ids))
    rows = cur.fetchall() or []
    dataset_ids = {r[0] for r in rows}
    cur.close()
    return dataset_ids

def get_fk_tables_referencing_dataset(conn: Connection) -> list:
    """
    Return list of (table_name, column_name) that have a foreign key referencing dataset(id).
    Uses INFORMATION_SCHEMA.KEY_COLUMN_USAGE.
    """

    cur = conn.get_cursor()
    sql = """
        SELECT TABLE_NAME, COLUMN_NAME
        FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE
        WHERE REFERENCED_TABLE_NAME = 'dataset'
          AND REFERENCED_COLUMN_NAME = 'id'
          AND CONSTRAINT_SCHEMA = %s
    """
    cur.execute(sql, ('gear_portal', ))
    rows = cur.fetchall() or []
    # mysql cursor default returns tuples; adapt if dictionary cursor used
    tables = []
    for r in rows:
        if isinstance(r, dict):
            tables.append((r["TABLE_NAME"], r["COLUMN_NAME"]))
        else:
            tables.append((r[0], r[1]))
    cur.close()
    # always include the dataset table itself
    if ("dataset", "id") not in tables:
        tables.insert(0, ("dataset", "id"))
    return tables

def get_fk_tables_not_referencing_dataset(conn: Connection) -> list:
    """
    Return list of (table_name, column_name) that have a foreign key the do not
    referencing dataset(id). Uses INFORMATION_SCHEMA.KEY_COLUMN_USAGE.
    """

    cur = conn.get_cursor()
    sql = """
        SELECT TABLE_NAME
        FROM INFORMATION_SCHEMA.TABLES
        WHERE TABLE_SCHEMA = %s
          AND TABLE_TYPE = 'BASE TABLE'
          AND TABLE_NAME NOT IN (
              SELECT DISTINCT TABLE_NAME
              FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE
              WHERE REFERENCED_TABLE_NAME = 'dataset'
                AND REFERENCED_COLUMN_NAME = 'id'
                AND CONSTRAINT_SCHEMA = %s
          )
    """
    cur.execute(sql, ("gear_portal", "gear_portal"))
    rows = cur.fetchall() or []
    # mysql cursor default returns tuples; adapt if dictionary cursor used
    tables = []
    for r in rows:
        if isinstance(r, dict):
            tables.append((r["TABLE_NAME"]))
        else:
            tables.append((r[0]))
    cur.close()
    return tables

def run_mysqldump_schema(servercfg: "configparser.ConfigParser", tables: list, out_path: str) -> bool:
    if not shutil.which("mysqldump"):
        print("mysqldump not found on PATH. Please run mysqldump manually.", file=sys.stderr)
        return False
    mysql_cfg = servercfg["database"] if servercfg and "database" in servercfg else {}
    user = mysql_cfg.get("user", "")
    passwd = mysql_cfg.get("password", "")
    host = mysql_cfg.get("host", "localhost")
    db = mysql_cfg.get("name", "gear_portal")

    cmd = ["mysqldump", "-h", host, "-u", user, "--no-data", db] + tables
    if passwd:
        cmd.insert(3, f"--password={passwd}")
    print("Dumping schema for tables:", tables)
    print(f"Running command: {' '.join(cmd)}")
    with open(out_path, "wb") as fh:
        proc = subprocess.run(cmd, stdout=fh)
    return proc.returncode == 0


def append_table_data_for_ids(servercfg: "configparser.ConfigParser", table: str, col: str, ids: list, out_path: str) -> int:
    """
    Append data rows for given ids for (table,col) to out_path.
    Returns number of rows written as reported by mysqldump returncode (not precise) -- we return len(ids) as best-effort.
    """
    if not ids:
        return 0
    if not shutil.which("mysqldump"):
        raise RuntimeError("mysqldump not found on PATH")

    mysql_cfg = servercfg["database"] if servercfg and "database" in servercfg else {}
    user = mysql_cfg.get("user", "")
    passwd = mysql_cfg.get("password", "")
    host = mysql_cfg.get("host", "localhost")
    db = mysql_cfg.get("name", "gear_portal")

    total_written = 0
    for chunk in chunked_iter(iter(ids), size=CHUNK):
        in_list = ",".join(sql_quote(x) for x in chunk)
        where_clause = f"{col} IN ({in_list})"
        cmd = [
            "mysqldump",
            "-h",
            host,
            "-u",
            user,
            "--single-transaction",
            "--quick",
            "--lock-tables=false",
            "--skip-add-drop-table",
            "--no-create-info",
            db,
            table,
            "--where",
            where_clause,
        ]
        if passwd:
            cmd.insert(3, f"--password={passwd}")
        print(f"Dumping data for {table} WHERE {col} IN (..{len(chunk)} ids..)")
        with open(out_path, "ab") as fh:
            proc = subprocess.run(cmd, stdout=fh)
        if proc.returncode != 0:
            print(f"Warning: mysqldump returned {proc.returncode} for table {table}", file=sys.stderr)
        total_written += len(chunk)
    return total_written


def main(argv: list | None = None):
    parser = argparse.ArgumentParser(description="Dump only dataset-related rows referenced by layouts.")
    parser.add_argument("--layout-ids", "-l", nargs="+", type=int, required=True, help="Layout id(s) to include")
    parser.add_argument("--dump-file", "-o", required=True, help="Output mysqldump filename")
    args = parser.parse_args(argv)

    servercfg = ServerConfig().parse()
    conn = Connection()

    TODO: Drop layouts that do not have any of the display IDs for these datasets

    try:
        dataset_ids = get_display_and_dataset_ids(conn, args.layout_ids)
        if not dataset_ids:
            print("No datasets found for given layout ids. Exiting.")
            return

        print(f"Found {len(dataset_ids)} dataset ids. List: {list(dataset_ids)[:10]}")

        fk_tables = get_fk_tables_referencing_dataset(conn)
        # fk_tables: list of (table, column)
        tables = [t for t, _ in fk_tables]
        # unique preserve order
        seen = set()
        tables_unique = [x for x in tables if not (x in seen or seen.add(x))]

        # Add all remaining tables to the dump.  We only want to restrict those tables that use datasets directly or indirectly
        not_fk_tables = get_fk_tables_not_referencing_dataset(conn)
        for t in not_fk_tables:
            if t not in seen:
                tables_unique.append(t)
                seen.add(t)

        # Step 1: dump schema for these tables
        success = run_mysqldump_schema(servercfg, tables_unique, args.dump_file)
        if not success:
            print("Schema dump failed. Aborting.", file=sys.stderr)
            return

        # Step 2: append data per table
        for table, col in fk_tables:
            # for dataset table col will be 'id'
            written = append_table_data_for_ids(servercfg, table, col, list(dataset_ids), args.dump_file)
            print(f"Appended ~{written} rows for {table}")

        for table in not_fk_tables:
            # dump all rows (col "1" in (1) always true)
            written = append_table_data_for_ids(servercfg, table, "1", [1], args.dump_file)  # where 1=1
            print(f"Appended ~{written} rows for {table}")

        print(f"Completed dump to {args.dump_file}")

    finally:
        try:
            conn.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()