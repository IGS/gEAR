from flask import request
from flask_restful import Resource, reqparse
import requests

DEFAULT_START_PAGE = 1
DEFAULT_PAGINATION = 20
DEFAULT_SORT = "file_id:asc"

# These are the facets updated on the UI
DEFAULT_FACETS = "sample.organism,sample.brain_sub_region_fullnames,sample.subspecimen_type,sample.technique,sample.modality"

# TODO: Add LOOM later
FILE_FORMATS_TO_INCLUDE = ["MEX", "H5AD", "H5"]
QUERIES_TO_INCLUDE = {
                        "file.format":FILE_FORMATS_TO_INCLUDE
                        , "access": "open"
                    }

PORTAL_FILES_API = "https://portal.nemoarchive.org/api/files"

def build_filters(query):
    """Build NeMO Archive portal API "filters" param based on supplied query

    Filters look like this
    {"op":"and","content":
        [
            {"op":"in","content":
                {"field":"file.format","value":
                    ["MEX","LOOM","H5"]
                }
            },{"op":"in","content":
                {"field":"subject.grant","value":
                    ["U01_Snyder-Mackler"]
                }
            }
        ]
    }

    Args:
        query (dict): Information used to filter datasets
    """

    filters = {"op": "and", "content": []}
    for field in query:
        field_filter = {"field": field, "value": query[field]}
        top_filter = {{"op": "in", "content": {field_filter}}}
        filters["content"].append(top_filter)

    return filters

def merge_query_and_defaults(query, default_queries):
    # merge both dicts, but keep "query" entries if fields in both dicts
    # (https://stackoverflow.com/a/26853961)
    merged_queries = default_queries | query

    # Now update fields that are in query and defaults
    for field in query:
        if field in QUERIES_TO_INCLUDE:
            if type(query[field]) == list:
                field_set = set(query[field])
                field_set.update(QUERIES_TO_INCLUDE[field])
                merged_queries[field] = list(field_set)
            else:
                pass    # Implement this later

    return merged_queries

class Query(Resource):
    def post(self):
        session_id = request.cookies.get('gear_session_id')
        parser = reqparse.RequestParser(bundle_errors=True)
        parser.add_argument('query', help='query filters dict required', type=dict, required=True)
        args = parser.parse_args()

        payload = {
            "facets": DEFAULT_FACETS
            , "from" : DEFAULT_START_PAGE
            , "size" : DEFAULT_PAGINATION
            , "sort" : DEFAULT_SORT
        }

        # Build the filters JSON that is used to query the portal API
        query = merge_query_and_defaults(args["query"], QUERIES_TO_INCLUDE)
        payload["filters"] = build_filters(query)

        # POST to NeMO Archive portal API

        # WARNING: Disabling SSL verification in the POST call
        # Throw error if things went awry (check apache ssl_error logs)
        result = requests.post(url, json=payload, verify=False)
        result.raise_for_status()

        return result.json()

"""
sample return value

{
    "data": {
        "hits": [
            {
                "file": {
                    "identifier": "nemo:der-heokos8",
                    "study": "MACOSKO_REGEV",
                    "mex_genes": "gs://nemo-public/biccn-unbundled/grant/u19_huang/macosko_regev/transcriptome/sncell/10X_v3/mouse/processed/counts/CellRanger5/pBICCNsMMrCBAN2iM002d190312b/pBICCNsMMrCBAN2iM002d190312b.features.tsv.gz",
                    "file_name": "pBICCNsMMrCBAN2iM002d190312b.mex.tar.gz",
                    "format": "MEX",
                    "mtime": 1619804353,
                    "mex_matrix": "gs://nemo-public/biccn-unbundled/grant/u19_huang/macosko_regev/transcriptome/sncell/10X_v3/mouse/processed/counts/CellRanger5/pBICCNsMMrCBAN2iM002d190312b/pBICCNsMMrCBAN2iM002d190312b.matrix.mtx.gz",
                    "node_type": "transcriptome",
                    "file": "https://data.nemoarchive.org/biccn/grant/u19_huang/macosko_regev/transcriptome/sncell/10X_v3/mouse/processed/counts/CellRanger5/pBICCNsMMrCBAN2iM002d190312b.mex.tar.gz",
                    "size": 37856007,
                    "subtype": "counts",
                    "id": "7b53d7eb-ca5a-450b-9f47-35c5574999a3",
                    "mex_barcodes": "gs://nemo-public/biccn-unbundled/grant/u19_huang/macosko_regev/transcriptome/sncell/10X_v3/mouse/processed/counts/CellRanger5/pBICCNsMMrCBAN2iM002d190312b/pBICCNsMMrCBAN2iM002d190312b.barcodes.tsv.gz",
                    "md5": "338dfcf9e7fb978c9a5cd865744e6d32",
                    "access": "open"
                },
                "file_id": "7b53d7eb-ca5a-450b-9f47-35c5574999a3"
            },
            ...
        ],
        "aggregations": {
            "file.format": {
                "buckets": [
                    {
                        "key": "FASTQ",
                        "doc_count": 528360
                    },
                    {
                        "key": "MEX",
                        "doc_count": 1417
                    },
                    {
                        "key": "BAM",
                        "doc_count": 399221
                    },
                    {
                        "key": "LIST",
                        "doc_count": 90
                    },
                    ...
                ]
            },
            "file.node_type": {
                "buckets": [
                    {
                        "key": "transcriptome",
                        "doc_count": 564
                    },
                    {
                        "key": "transcriptomics",
                        "doc_count": 781
                    },
                    {
                        "key": "multimodal",
                        "doc_count": 4
                    },
                    {
                        "key": "epigenome",
                        "doc_count": 110
                    }
                ]
            },
            "file.subtype": {
                "buckets": [
                    {
                        "key": "counts",
                        "doc_count": 1449
                    },
                    {
                        "key": "unknown",
                        "doc_count": 10
                    }
                ]
            }
        },
        "pagination": {
            "count": 20,
            "total": 1459,
            "page": 1,
            "pages": 73,
            "from": 1,
            "sort": "file_id:asc",
            "size": 20,
            "sample_total": 1428
        }
    }
}
"""