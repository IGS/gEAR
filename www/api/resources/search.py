"""
search.py - The unified search across gEAR entities (#1383, Stage 1).

Serves /search in www/api/api.py. The searching itself is in lib/gear/search.

GET parameters (all optional): q, types, organism_ids, ownership, date_added, sort, page, limit
(see gear.search.parse_criteria). The user is identified by the gear_session_id cookie.

Responses:
- 200: {"success": 1, "problem": "", "counts": {...}, "results": [...], "pagination": {...}}
- 400: {"success": 0, "problem": "<what was wrong with the request>"}
- 500: {"success": 0, "problem": "<a general message>"}; the details go to the server log only
"""

import sys
import traceback

import geardb
from flask import request
from flask_restful import Resource
from gear.search import parse_criteria, run_search
from gear.search.cell_cell import CellCellSearcher
from gear.search.datasets import DatasetSearcher
from gear.search.gene_lists import GeneListSearcher
from gear.search.layouts import CollectionSearcher
from gear.search.projections import ProjectionSearcher
from gear.search.study_pages import StudyPageSearcher

SEARCH_FAILED = "The search could not be run. Please try different search terms."


def build_searchers():
    """One searcher per result type in gear.search.TYPES."""
    return {
        "dataset": DatasetSearcher(),
        "collection": CollectionSearcher(),
        "genelist": GeneListSearcher(),
        # Only this site's study pages are searched
        "study": StudyPageSearcher(geardb.domain_label),
        "projection": ProjectionSearcher(),
        "ccc": CellCellSearcher(),
    }


class UnifiedSearch(Resource):
    def get(self):
        try:
            criteria = parse_criteria(request.args)
        except ValueError as err:
            return {"success": 0, "problem": str(err)}, 400

        session_id = request.cookies.get("gear_session_id")
        user = geardb.get_user_from_session_id(session_id) if session_id else None

        cnx = geardb.Connection()
        cursor = cnx.get_cursor()
        try:
            result = run_search(criteria, user, cursor, build_searchers())
        except Exception:
            print(f"search.py: search failed for {dict(request.args)}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            return {"success": 0, "problem": SEARCH_FAILED}, 500
        finally:
            cursor.close()
            cnx.close()

        return {"success": 1, "problem": "", **result}
