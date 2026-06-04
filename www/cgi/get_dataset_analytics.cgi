#!/opt/bin/python3

"""
Returns comprehensive database analytics including:
- Total database size
- Number of datasets (by type, public/private, status)
- User and organism statistics
- Dataset size info
"""

import cgi
import json
import os
import sys
from pathlib import Path
from datetime import datetime

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

import geardb


def get_directory_size(path):
    """Calculate total size of directory in GB"""
    if not os.path.exists(path):
        return 0
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            try:
                total_size += os.path.getsize(filepath)
            except (OSError, FileNotFoundError):
                pass
    return total_size / (1024 ** 3)  # Convert to GB


def get_h5ad_gene_count(h5ad_path):
    """Get number of genes from H5AD file"""
    if not os.path.exists(h5ad_path):
        return 0

    try:
        import h5py
        with h5py.File(h5ad_path, 'r') as f:
            # H5AD stores genes in 'var' (variables)
            if 'var' in f:
                return len(f['var'])
            elif 'X' in f:
                # Alternative: get from data matrix shape
                return f['X'].shape[1]
    except Exception as e:
        print(f"Error reading H5AD file {h5ad_path}: {e}", file=sys.stderr)

    return 0


def main():
    try:
        cnx = geardb.Connection()
        cursor = cnx.get_cursor()

        # Get total number of datasets
        cursor.execute("""
            SELECT COUNT(*) as total_count
            FROM dataset
            WHERE marked_for_removal = 0
        """)
        result = cursor.fetchone()
        total_datasets = result[0] if result else 0

        # Get public vs private dataset counts
        cursor.execute("""
            SELECT is_public, COUNT(*) as count
            FROM dataset
            WHERE marked_for_removal = 0
            GROUP BY is_public
        """)
        public_datasets = 0
        private_datasets = 0
        for row in cursor.fetchall():
            if row[0] == 1:
                public_datasets = row[1]
            else:
                private_datasets = row[1]

        # Get dataset count by type
        cursor.execute("""
            SELECT dtype, COUNT(*) as count
            FROM dataset
            WHERE marked_for_removal = 0
            GROUP BY dtype
        """)

        dtype_counts = {}
        dataset_types = {}

        for row in cursor.fetchall():
            dtype = row[0] if row[0] else 'unknown'
            count = row[1]
            dtype_counts[dtype] = count
            dataset_types[dtype] = {'count': count}

        # Get dataset status breakdown
        cursor.execute("""
            SELECT load_status, COUNT(*) as count
            FROM dataset
            WHERE marked_for_removal = 0
            GROUP BY load_status
        """)
        status_breakdown = {}
        for row in cursor.fetchall():
            status = row[0] if row[0] else 'unknown'
            status_breakdown[status] = row[1]

        # Get number of unique organisms
        cursor.execute("""
            SELECT COUNT(DISTINCT organism_id) as organism_count
            FROM dataset
            WHERE marked_for_removal = 0
        """)
        result = cursor.fetchone()
        organism_count = result[0] if result else 0

        # Get number of registered users
        cursor.execute("""
            SELECT COUNT(*) as user_count
            FROM guser
        """)
        result = cursor.fetchone()
        user_count = result[0] if result else 0

        # Get number of shared datasets
        cursor.execute("""
            SELECT COUNT(DISTINCT dataset_id) as shared_count
            FROM dataset_shares
        """)
        result = cursor.fetchone()
        shared_dataset_count = result[0] if result else 0

        # Get all dates when datasets were added
        cursor.execute("""
            SELECT DATE(date_added) as date, COUNT(*) as count
            FROM dataset
            WHERE marked_for_removal = 0 AND date_added IS NOT NULL
            GROUP BY DATE(date_added)
            ORDER BY DATE(date_added) DESC
        """)
        dates_added = {}
        for row in cursor.fetchall():
            date = row[0].isoformat() if row[0] else None
            if date:
                dates_added[date] = row[1]

        # Get number of gene lists
        cursor.execute("""
            SELECT COUNT(*) as gene_cart_count
            FROM gene_cart
        """)
        result = cursor.fetchone()
        gene_cart_count = result[0] if result else 0

        # Get number of layouts
        cursor.execute("""
            SELECT COUNT(*) as layout_count
            FROM layout
        """)
        result = cursor.fetchone()
        layout_count = result[0] if result else 0

        # Get total number of dataset shares
        cursor.execute("""
            SELECT COUNT(*) as total_shares
            FROM dataset_shares
        """)
        result = cursor.fetchone()
        total_shares = result[0] if result else 0

        # Calculate database size
        total_size_gb = 0

        # Check multiple dataset directories
        dataset_base_paths = [
            os.path.abspath(os.path.join('..', 'datasets')),
            os.path.abspath(os.path.join('..', 'datasets_uploaded')),
            os.path.abspath(os.path.join('..', 'datasets_epigenetic'))
        ]

        # Size of all uploaded datasets
        for datasets_dir in dataset_base_paths:
            if os.path.exists(datasets_dir):
                try:
                    total_size_gb += get_directory_size(datasets_dir)
                except Exception as e:
                    print(f"Error calculating size for {datasets_dir}: {e}", file=sys.stderr)

        # Size of database (MySQL data directory estimate)
        # Try multiple common locations
        mysql_data_dirs = [
            "/var/lib/mysql",
            "/usr/local/var/mysql",
            os.path.expanduser("~/mysql"),
            "/opt/mysql"
        ]

        for mysql_data_dir in mysql_data_dirs:
            if os.path.exists(mysql_data_dir):
                # Try to get size of gear database specifically
                gear_db_dir = os.path.join(mysql_data_dir, "gear")
                if os.path.exists(gear_db_dir):
                    try:
                        total_size_gb += get_directory_size(gear_db_dir)
                        break
                    except Exception as e:
                        print(f"Error reading MySQL directory {gear_db_dir}: {e}", file=sys.stderr)
                        continue

        # Calculate average dataset size
        avg_dataset_size = round(total_size_gb / max(total_datasets, 1), 4)

        # Prepare response
        response = {
            'overview': {
                'totalSize': round(total_size_gb, 2),
                'numberOfDatasets': total_datasets,
                'averageDatasetSize': avg_dataset_size,
                'publicDatasets': public_datasets,
                'privateDatasets': private_datasets
            },
            'byType': {},
            'statusBreakdown': status_breakdown,
            'statistics': {
                'totalOrganisms': organism_count,
                'registeredUsers': user_count,
                'sharedDatasets': shared_dataset_count,
                'totalShares': total_shares,
                'geneLists': gene_cart_count,
                'layouts': layout_count,
                'datesAdded': dates_added
            }
        }

        # Map dtype values to friendly names
        type_mapping = {
            'svg-expression': 'Microarray',
            'h5ad-bulk': 'RNA-Seq',
            'h5ad-singlecell': 'Single-Cell',
            'h5ad-epigenetic': 'Epigenetic',
            'unknown': 'Other'
        }

        for dtype, data in dataset_types.items():
            friendly_name = type_mapping.get(dtype, dtype)
            response['byType'][friendly_name] = {
                'count': data['count']
            }

        cursor.close()
        cnx.close()

        print('Content-Type: application/json\n\n')
        print(json.dumps(response))

    except Exception as e:
        print('Content-Type: application/json\n\n')
        print(json.dumps({
            'error': str(e),
            'overview': {
                'totalSize': 0,
                'numberOfDatasets': 0,
                'averageDatasetSize': 0,
                'publicDatasets': 0,
                'privateDatasets': 0
            },
            'byType': {},
            'statusBreakdown': {},
            'statistics': {}
        }))


if __name__ == '__main__':
    main()
