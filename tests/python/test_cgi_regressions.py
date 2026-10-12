"""Regression tests for closed GitHub issues in CGI scripts (one class per issue or area)."""

import re
import shutil
import uuid
from pathlib import Path

import pytest
from helpers.cgi_harness import CGI_DIR, run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
USERHISTORY = {"gear.userhistory": FAKES / "accounts_userhistory.py"}
SESSIONS = {"owner": 1, "other": 2}


def assert_clean_exit(result):
    assert result.returncode == 0, result.stderr[-2000:]
    assert "Traceback" not in result.stderr, result.stderr[-2000:]


# ---------------------------------------------------------------- get_shared_info.cgi

class TestGetSharedInfo:
    EMPTY = {"gene_symbol": None, "dataset_label": None, "layout_label": None, "owner_name": None}
    DB = {
        "users": {"1": {"user_name": "Sally Nguyen"}},
        "datasets": [{"id": "DS1", "share_id": "S1", "owner_id": 1, "title": "Utricle",
                      "owner_display": {"plotly_config": '{"gene_symbol": "Sox2"}'}}],
        "layouts": [{"id": 10, "share_id": "L1", "user_id": 1, "label": "My collection",
                     "members": [{"dataset_id": "DS1"}]},
                    {"id": 11, "share_id": "L2", "user_id": 1, "label": "Empty collection"}],
    }

    @pytest.mark.parametrize("query", [{"dataset_share_id": "NOPE"}, {"layout_share_id": "NOPE"},
                                       {"dataset_share_id": "null", "layout_share_id": "null"}, {}])
    def test_unknown_share_id_returns_empty_info(self, query):
        """Regression for #1001: an unknown dataset or collection share ID returns empty JSON, not an UnboundLocalError."""
        result = run_cgi("get_shared_info.cgi", db=self.DB, query=query)
        assert_clean_exit(result)
        assert result.json() == self.EMPTY

    def test_dataset_share(self):
        """Regression for #1001: a known dataset share returns its title, owner and the owner's curated gene."""
        result = run_cgi("get_shared_info.cgi", db=self.DB, query={"dataset_share_id": "S1"})
        assert result.json() == {**self.EMPTY, "dataset_label": "Utricle", "owner_name": "Sally Nguyen",
                                 "gene_symbol": "Sox2"}

    def test_layout_share_uses_first_dataset(self):
        """Regression for #1001: a known collection share returns its label, owner and its first dataset's gene."""
        result = run_cgi("get_shared_info.cgi", db=self.DB, query={"layout_share_id": "L1"})
        assert result.json() == {**self.EMPTY, "layout_label": "My collection", "owner_name": "Sally Nguyen",
                                 "gene_symbol": "Sox2"}

    def test_empty_layout_share(self):
        """Regression for #1001: a collection with no datasets still returns its label and owner."""
        result = run_cgi("get_shared_info.cgi", db=self.DB, query={"layout_share_id": "L2"})
        assert_clean_exit(result)
        assert result.json() == {**self.EMPTY, "layout_label": "Empty collection", "owner_name": "Sally Nguyen"}


# ---------------------------------------------------------------- search_genes.cgi

class TestSearchGenesHistory:
    DB = {"sessions": SESSIONS, "genes": [{"gene_symbol": "Sox2", "organism_id": 1}],
          "layouts": [{"id": 10, "share_id": "L1", "user_id": 1, "label": "Hearing"}]}

    def run(self, **query):
        return run_cgi("search_genes.cgi", db=self.DB, fake_modules=USERHISTORY,
                       query={"search_gene_symbol": "Sox2", "exact_match": "true", **query})

    @pytest.mark.parametrize("layout_share_id", ["NOPE", None])
    def test_unknown_layout_writes_no_history(self, layout_share_id):
        """Regression for #940: a logged-in search with an unknown layout_share_id logs no history and doesn't crash."""
        query = {"session_id": "owner"}
        if layout_share_id:
            query["layout_share_id"] = layout_share_id
        result = self.run(**query)
        assert_clean_exit(result)
        assert result.json() == {"sox2": {"by_organism": {"1": [{"gene_symbol": "Sox2", "organism_id": 1}]}}}
        assert not result.logged("userhistory.add_record")

    def test_known_layout_writes_history(self):
        """Regression for #940: a logged-in search with a known collection is still recorded in the user's history."""
        result = self.run(session_id="owner", layout_share_id="L1")
        assert_clean_exit(result)
        (record,) = result.logged("userhistory.add_record")
        assert record["user_id"] == 1 and record["layout_share_id"] == "L1"
        assert record["entry_category"] == "gene_search" and record["label"] == '"Sox2" in Hearing'

    def test_anonymous_search_writes_no_history(self):
        """Regression for #940: an anonymous search logs no history."""
        result = self.run(layout_share_id="L1")
        assert_clean_exit(result)
        assert "sox2" in result.json()
        assert not result.logged("userhistory.add_record")


# ---------------------------------------------------------------- get_gene_cart_members.cgi

@pytest.mark.parametrize("session", [None, "expired"])
def test_gene_cart_members_without_login(session):
    """Regression for #963: an anonymous (or unrecognised-session) request for a public gene list returns its genes."""
    db = {"sessions": SESSIONS,
          "gene_carts": [{"id": 5, "share_id": "GC1", "user_id": 1, "gctype": "unweighted-list", "is_public": 1}],
          "sql": [{"match": "FROM gene_cart_member", "rows": [[11, "Sox2"]]}]}
    query = {"share_id": "GC1"}
    if session:
        query["session_id"] = session
    result = run_cgi("get_gene_cart_members.cgi", db=db, query=query)
    assert result.json() == {"gene_symbols": [{"id": 11, "label": "Sox2"}], "success": 1}
    assert result.writes() == []


# ---------------------------------------------------------------- update_share_id.cgi (gene lists)

class TestUpdateGeneCartShareId:
    """Runs a copy of update_share_id.cgi from a temporary www/ so cart and projection files are disposable."""

    OLD, NEW = "ab12cd34", "my-markers"

    @pytest.fixture
    def www(self, tmp_path):
        www = tmp_path / "www"
        (www / "cgi").mkdir(parents=True)
        shutil.copy2(CGI_DIR / "update_share_id.cgi", www / "cgi" / "update_share_id.cgi")
        for sub in ("carts", "projections/by_dataset", "projections/by_genecart"):
            (www / sub).mkdir(parents=True)
        (www / "carts" / f"cart.{self.OLD}.tab").write_text("id\tgene\tP1\n")
        (www / "carts" / f"cart.{self.OLD}.h5ad").write_bytes(b"h5ad")
        return www

    def run(self, www, session_id="owner", gctype="weighted-list"):
        db = {"sessions": SESSIONS,
              "gene_carts": [{"id": 6, "share_id": self.OLD, "user_id": 1, "gctype": gctype}]}
        return run_cgi(str(www / "cgi" / "update_share_id.cgi"), db=db, query={
            "session_id": session_id, "scope": "genecart", "share_id": self.OLD, "new_share_id": self.NEW})

    def cart_files(self, www):
        return sorted(p.name for p in (www / "carts").iterdir())

    def test_weighted_cart_files_are_renamed(self, www):
        """Regression for #978: renaming a weighted gene list's share ID also renames its cart.<id>.tab/.h5ad files."""
        result = self.run(www)
        assert_clean_exit(result)
        assert result.json() == {"error": "", "success": 1}
        assert result.logged("gene_cart.save_change", self.NEW)
        assert self.cart_files(www) == [f"cart.{self.NEW}.h5ad", f"cart.{self.NEW}.tab"]

    def test_other_carts_containing_the_old_id_are_not_renamed(self, www):
        """Regression for #978: only this list's own files are renamed, not other carts whose names contain its ID."""
        (www / "carts" / f"cart.{self.OLD}99.tab").write_text("other list\n")
        (www / "carts" / f"cart.x{self.OLD}.h5ad").write_bytes(b"other list")
        assert self.run(www).json()["success"] == 1
        assert self.cart_files(www) == sorted([f"cart.{self.NEW}.h5ad", f"cart.{self.NEW}.tab",
                                               f"cart.{self.OLD}99.tab", f"cart.x{self.OLD}.h5ad"])

    def test_non_owner_is_refused(self, www):
        """Regression for #978: another user can't rename the list or its files."""
        result = self.run(www, session_id="other")
        assert result.json() == {"error": "User does not own this genecart.", "success": 0}
        assert not result.logged("gene_cart.save_change")
        assert result.writes() == []
        assert self.cart_files(www) == [f"cart.{self.OLD}.h5ad", f"cart.{self.OLD}.tab"]


# ---------------------------------------------------------------- remove_dataset.cgi

class TestRemoveDataset:
    DB = {"sessions": SESSIONS, "sql": [{"match": "FROM dataset d WHERE d.id", "rows": [["DS1", 1]]}]}

    def test_owner_removes_dataset(self):
        """Regression for #860: deleting from layout_displays uses multi-table DELETE syntax (it joins dataset_display)."""
        result = run_cgi("remove_dataset.cgi", db=self.DB, query={"session_id": "owner", "dataset_id": "DS1"})
        assert result.json() == {"success": 1, "dataset": []}
        layout_displays, shares, dataset = result.writes()
        assert layout_displays["query"].startswith(
            "DELETE layout_displays FROM layout_displays JOIN dataset_display "
            "ON layout_displays.display_id = dataset_display.id")
        assert layout_displays["params"] == ["DS1"]
        assert shares["query"] == "DELETE FROM dataset_shares WHERE dataset_id = %s" and shares["params"] == ["DS1"]
        assert dataset["query"].startswith("UPDATE dataset SET marked_for_removal = 1") and dataset["params"] == ["DS1"]
        assert result.logged("commit")

    @pytest.mark.parametrize("session, dataset_id", [("other", "DS1"), (None, "DS1"), ("owner", "NOPE")])
    def test_non_owner_is_refused(self, session, dataset_id):
        """Regression for #860: someone who doesn't own the dataset (or an unknown dataset) gets an error and no writes."""
        db = self.DB if dataset_id == "DS1" else {"sessions": SESSIONS}
        query = {"dataset_id": dataset_id}
        if session:
            query["session_id"] = session
        result = run_cgi("remove_dataset.cgi", db=db, query=query)
        assert result.json() == {"success": 0,
                                 "error": "Not able to remove dataset. User does not own the dataset."}
        assert result.writes() == []
        assert not result.logged("commit")


# ---------------------------------------------------------------- delete_dataset_analysis.cgi

class TestDeleteDatasetAnalysis:
    """Uses the real gear.analysis with its root moved to tmp_path (fakes/regress_analysis.py)."""

    DATASET_ID = "5f0c1a6e-0000-4000-8000-000000000001"
    ANALYSIS_ID = "9d7e2b4c-1111-4111-8111-111111111111"
    DB = {"sessions": SESSIONS, "datasets": [{"id": DATASET_ID, "owner_id": 1}]}

    def run(self, root, session_id="owner", **form):
        form = {"session_id": session_id, "dataset_id": self.DATASET_ID, **form}
        return run_cgi("delete_dataset_analysis.cgi", db=self.DB, form=form, env={"REGRESS_ROOT_DIR": str(root)},
                       fake_modules={"gear.analysis": FAKES / "regress_analysis.py"})

    def saved_dir(self, root, user_id, analysis_id):
        path = root / "www" / "analyses" / "by_user" / str(user_id) / self.DATASET_ID / analysis_id
        path.mkdir(parents=True)
        (path / f"{self.DATASET_ID}.h5ad").write_bytes(b"h5ad")
        return path

    @pytest.fixture
    def primary(self, tmp_path):
        datasets = tmp_path / "www" / "datasets"
        datasets.mkdir(parents=True)
        (datasets / f"{self.DATASET_ID}.h5ad").write_bytes(b"primary")
        return datasets

    def test_owner_deletes_saved_analysis(self, tmp_path):
        """Regression for #1151: a user's saved analysis (UUID-named directory) is deleted, and nothing else."""
        target = self.saved_dir(tmp_path, 1, self.ANALYSIS_ID)
        sibling = self.saved_dir(tmp_path, 1, str(uuid.uuid4()))
        result = self.run(tmp_path, analysis_id=self.ANALYSIS_ID, analysis_type="user_saved")
        assert result.json() == {"success": 1}
        assert not target.exists() and sibling.exists()

    def test_missing_type_cannot_delete_primary_files(self, tmp_path, primary):
        """Regression for #1151: with no analysis_type (so a primary analysis), the www/datasets directory is not deleted."""
        result = self.run(tmp_path, analysis_id=self.ANALYSIS_ID)
        assert result.json()["success"] == 0
        assert (primary / f"{self.DATASET_ID}.h5ad").read_bytes() == b"primary"

    def test_primary_type_is_refused(self, tmp_path, primary):
        """Regression for #1151: an explicit primary analysis is refused."""
        result = self.run(tmp_path, analysis_id=self.DATASET_ID, analysis_type="primary")
        assert result.json() == {"success": 0, "error": "Cannot delete primary analyses"}
        assert (primary / f"{self.DATASET_ID}.h5ad").exists()

    @pytest.mark.parametrize("analysis_id", ["my-analysis", f"{ANALYSIS_ID}x"])
    def test_non_uuid_directory_is_refused(self, tmp_path, analysis_id):
        """Regression for #1151: a directory whose name is not a UUID is never deleted."""
        target = self.saved_dir(tmp_path, 1, analysis_id)
        result = self.run(tmp_path, analysis_id=analysis_id, analysis_type="user_saved")
        assert result.json()["success"] == 0
        assert target.exists()

    def test_path_in_analysis_id_cannot_reach_another_users_analysis(self, tmp_path):
        """Regression for #1151: analysis_id "../../2/<dataset>/<uuid>" can't delete user 2's analysis (its basename is a UUID)."""
        self.saved_dir(tmp_path, 1, str(uuid.uuid4()))    # user 1's own dataset directory exists
        victim = self.saved_dir(tmp_path, 2, self.ANALYSIS_ID)
        result = self.run(tmp_path, analysis_id=f"../../2/{self.DATASET_ID}/{self.ANALYSIS_ID}",
                          analysis_type="user_saved")
        assert result.json()["success"] == 0
        assert victim.exists()

    def test_public_analysis_cannot_be_deleted(self, tmp_path):
        """Regression for #1151: only the user's own saved/unsaved analyses can be deleted, not a public one."""
        public = tmp_path / "www" / "analyses" / "by_dataset" / self.DATASET_ID / self.ANALYSIS_ID
        public.mkdir(parents=True)
        result = self.run(tmp_path, session_id="other", analysis_id=self.ANALYSIS_ID, analysis_type="public")
        assert result.json()["success"] == 0
        assert public.exists()

    def test_other_users_analysis_is_not_found(self, tmp_path):
        """Regression for #1151: another user's request resolves to their own (missing) directory, so nothing is deleted."""
        target = self.saved_dir(tmp_path, 1, self.ANALYSIS_ID)
        result = self.run(tmp_path, session_id="other", analysis_id=self.ANALYSIS_ID, analysis_type="user_saved")
        assert result.json() == {"success": 0, "error": "Analysis not found"}
        assert target.exists()

    def test_anonymous_is_refused(self, tmp_path):
        """Regression for #1151: an anonymous request deletes nothing."""
        target = self.saved_dir(tmp_path, 1, self.ANALYSIS_ID)
        result = self.run(tmp_path, session_id="", analysis_id=self.ANALYSIS_ID, analysis_type="user_saved")
        assert result.json() == {"success": 0, "error": "User not found"}
        assert target.exists()


# ---------------------------------------------------------------- get_h5ad_dataset_list.cgi

class TestH5adDatasetListShareId:
    DB = {"sessions": SESSIONS, "datasets": [
        {"id": "PUB", "share_id": "SPUB", "owner_id": 1, "is_public": 1},
        {"id": "PRIV", "share_id": "SPRIV", "owner_id": 1, "is_public": 0},
        {"id": "ODD", "share_id": "null", "owner_id": 1, "is_public": 0},
    ]}

    def ids(self, body):
        return {key: [d["id"] for d in body[key]["datasets"]] for key in ("user", "public", "shared_with_user")}

    @pytest.mark.parametrize("session", [None, "other"])
    def test_share_id_adds_a_dataset_the_user_cant_otherwise_see(self, session):
        """Regression for #1010/#1094: a share_id for a private dataset adds it to "shared_with_user"."""
        query = {"share_id": "SPRIV"}
        if session:
            query["session_id"] = session
        result = run_cgi("get_h5ad_dataset_list.cgi", db=self.DB, query=query)
        assert self.ids(result.json()) == {"user": [], "public": ["PUB"], "shared_with_user": ["PRIV"]}

    def test_share_id_already_listed_is_not_duplicated(self):
        """Regression for #1010/#1094: a share_id for a dataset already listed isn't added again."""
        result = run_cgi("get_h5ad_dataset_list.cgi", db=self.DB, query={"session_id": "owner", "share_id": "SPUB"})
        assert self.ids(result.json()) == {"user": ["PUB", "PRIV", "ODD"], "public": ["PUB"], "shared_with_user": []}

    @pytest.mark.parametrize("share_id", ["null", "NOPE"])
    def test_null_or_unknown_share_id_adds_nothing(self, share_id):
        """Regression for #1010/#1094: share_id "null" is treated as no share_id, and an unknown one adds nothing."""
        result = run_cgi("get_h5ad_dataset_list.cgi", db=self.DB, query={"session_id": "other", "share_id": share_id})
        assert self.ids(result.json()) == {"user": [], "public": ["PUB"], "shared_with_user": []}


# ---------------------------------------------------------------- search_datasets.cgi / search_gene_carts.cgi

GROUP_SUBQUERY = re.compile(
    r"IN \(SELECT DISTINCT user_id FROM user_group_membership WHERE group_id IN "
    r"\(SELECT group_id FROM user_group_membership WHERE user_id = %s\)\)")
EMPTY_PARENS = re.compile(r"\(\s*\)")


def main_query(result, table):
    (entry,) = [e for e in result.logged("sql") if f"FROM {table}" in e["query"] and "COUNT(" not in e["query"]]
    return entry


@pytest.mark.parametrize("script, table, key", [
    ("search_datasets.cgi", "dataset d", "datasets"),
    ("search_gene_carts.cgi", "gene_cart gc", "gene_carts"),
])
class TestSearchOwnership:
    def run(self, script, **query):
        return run_cgi(script, db={"sessions": SESSIONS}, fake_modules=USERHISTORY, query=query)

    def test_group_ownership_with_no_matches_is_empty(self, script, table, key):
        """Regression for #787: ownership=group adds the user_group_membership subquery and returns no results, not an error."""
        result = self.run(script, session_id="owner", ownership="group")
        assert_clean_exit(result)
        body = result.json()
        assert body["success"] == 1 and body[key] == [] and body["pagination"]["total_results"] == 0
        entry = main_query(result, table)
        assert GROUP_SUBQUERY.search(entry["query"])
        assert entry["params"] == [1]

    def test_shared_only_ownership_is_empty_not_an_error(self, script, table, key):
        """Regression for #787: ownership=shared alone (nothing shared) returns no results instead of an SQL error."""
        result = self.run(script, session_id="owner", ownership="shared")
        body = result.json()
        assert body["success"] == 1 and body[key] == []
        assert not EMPTY_PARENS.search(main_query(result, table)["query"])

    def test_anonymous_ownership_group_sees_only_public(self, script, table, key):
        """Regression for #787: ownership is ignored for anonymous users, who only see public items."""
        result = self.run(script, ownership="group")
        assert result.json()["success"] == 1
        entry = main_query(result, table)
        assert "user_group_membership" not in entry["query"] and "is_public = 1" in entry["query"]


class TestSearchDatasetsSortOrder:
    def run(self, **query):
        db = {"sql": [{"match": "SELECT d.id", "rows": [["B"], ["A"], ["C"]]},
                      {"match": "SELECT COUNT(*)", "rows": [[3]]}]}
        return run_cgi("search_datasets.cgi", db=db, fake_modules=USERHISTORY, query=query)

    def test_sort_by_owner_joins_guser(self):
        """Regression for #756/#781: sort_by=owner joins guser on the owner and orders by g.user_name, in both queries."""
        result = self.run(sort_by="owner")
        assert result.json()["pagination"]["total_results"] == 3
        main = main_query(result, "dataset d")["query"]
        assert "SELECT d.id, g.user_name FROM dataset d, guser g" in main
        assert "d.owner_id = g.id" in main and re.search(r"ORDER BY\s+g\.user_name", main)
        (count,) = result.logged("sql", "COUNT(*)")
        assert "FROM dataset d, guser g" in count["query"] and "d.owner_id = g.id" in count["query"]

    @pytest.mark.parametrize("sort_by", ["owner", "title", "date"])
    def test_result_order_is_kept(self, sort_by):
        """Regression for #756: dataset IDs are looked up in the order the sorted query returned them."""
        result = self.run(sort_by=sort_by)
        (lookup,) = result.logged("dataset_collection.get_by_dataset_ids")
        assert lookup["ids"] == ["B", "A", "C"]
