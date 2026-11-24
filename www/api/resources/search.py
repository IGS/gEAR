"""
search.py - The unified search across gEAR entities (#1383, Stage 1).

Serves /search in www/api/api.py. The searching itself is in lib/gear/search.

GET parameters (all optional): q, types, organism_ids, ownership, date_added, sort, page, limit
(see gear.search.parse_criteria). The user is identified by the gear_session_id cookie.
A logged-in user's new keyword searches are added to their history (see log_search).

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
from gear.userhistory import UserHistory

SEARCH_FAILED = "The search could not be run. Please try different search terms."
HISTORY_CATEGORY = "unified_search"


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


def log_search(user, criteria, cnx):
    """
    Add a logged-in user's keyword search to their history (the recent activity list on the home
    page), as the Dataset Explorer and Gene List Manager searches do.

    The page reruns the search on every filter, sort or page change, so only a new search is
    logged: page 1, with keywords that differ from the user's latest history entry. A failure
    here is reported in the server log and never fails the search.
    """
    if not user or not criteria.q or criteria.page > 1:
        return

    label = f"Results matching '{criteria.q}'"
    try:
        history = UserHistory(connection=cnx, user_id=user.id)
        if not history.user_history_enabled:
            return
        latest = history.get_latest_entries(entry_count=1)
        if latest and latest[0].entry_category == HISTORY_CATEGORY and latest[0].label == label:
            return
        history.add_record(user_id=user.id, entry_category=HISTORY_CATEGORY, label=label,
                           search_terms=criteria.q.split())
    except Exception:
        print("search.py: could not add the search to the user's history", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)


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
            log_search(user, criteria, cnx)
        except Exception:
            print(f"search.py: search failed for {dict(request.args)}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            return {"success": 0, "problem": SEARCH_FAILED}, 500
        finally:
            cursor.close()
            cnx.close()

        return {"success": 1, "problem": "", **result}
