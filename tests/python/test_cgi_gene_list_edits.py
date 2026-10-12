"""Gene list CGIs: creating lists (form and JSON), editing them, listing them and previewing their contents."""

import io
import json
import shutil
import subprocess
from pathlib import Path
from unittest import mock

import anndata
import numpy as np
import pandas as pd
import pytest
from helpers.cgi_harness import CGI_DIR, run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
USERHISTORY = {"gear.userhistory": FAKES / "cgi_acct_userhistory.py"}

# Session "owner" is user 1; "other" is user 2.
SESSIONS = {"owner": 1, "other": 2}

WEIGHTED_TAB = "id\tgene\tP1\tP2\nE1\tSox2\t0.5\t1\nE2\tPax6\t-1.5\t0\nE3\tAtoh1\t2.0\t1\n"


def assert_clean_exit(result):
    assert result.returncode == 0, result.stderr[-2000:]
    assert "Traceback" not in result.stderr, result.stderr[-2000:]


@pytest.fixture
def www(tmp_path):
    """A temporary www/ with an empty carts/ directory; copy_cgi() puts a CGI in its cgi/ directory."""
    www = tmp_path / "www"
    (www / "cgi").mkdir(parents=True)
    (www / "carts").mkdir()
    return www


def copy_cgi(www, script):
    shutil.copy2(CGI_DIR / script, www / "cgi" / script)
    return str(www / "cgi" / script)


def cart_files(www):
    return sorted(p.name for p in (www / "carts").iterdir())


def run_json_cgi(script, payload, **kwargs):
    """Run a CGI that reads a JSON request body from stdin (run_cgi only sends query strings and forms)."""
    real_run = subprocess.run

    def with_body(*args, **run_kwargs):
        return real_run(*args, **{**run_kwargs, "input": json.dumps(payload).encode()})

    with mock.patch("helpers.cgi_harness.subprocess.run", with_body):
        return run_cgi(script, env={"REQUEST_METHOD": "POST", "CONTENT_TYPE": "application/json"}, **kwargs)


# ---------------------------------------------------------------- save_new_genecart_form.cgi

@pytest.mark.slow   # weighted lists are converted with scanpy in the CGI subprocess
class TestSaveNewGeneCartForm:
    DB = {"sessions": SESSIONS, "lastrowid": 7}
    FIELDS = {"new_cart_label": "Hair cells", "new_cart_organism_id": "1", "new_cart_ldesc": "Markers",
              "is_public": "0"}

    def run(self, www, session_id="owner", files=None, **form):
        form = {**self.FIELDS, **form}
        if session_id:
            form["session_id"] = session_id
        return run_cgi(copy_cgi(www, "save_new_genecart_form.cgi"), db=self.DB, fake_modules=USERHISTORY,
                       form=form, files=files)

    def saved(self, result):
        (save,) = result.logged("gene_cart.save")
        return save

    def assert_refused(self, result, www, message):
        assert result.status == 500
        assert message in result.json()["message"]
        assert not result.logged("gene_cart.save")
        assert not result.logged("userhistory.add_record")
        assert cart_files(www) == []

    @pytest.mark.parametrize("session_id", [None, "expired"])
    def test_anonymous_is_refused(self, www, session_id):
        result = self.run(www, session_id=session_id, new_cart_upload_type="pasted_genes",
                          new_cart_pasted_genes="Sox2")
        self.assert_refused(result, www, "No logged-in user")

    def test_pasted_genes(self, www):
        result = self.run(www, new_cart_upload_type="pasted_genes", new_cart_pasted_genes="Sox2, Pax6 Sox2\nAtoh1")
        assert_clean_exit(result)
        assert result.json() == {"id": 7}
        save = self.saved(result)
        assert save["genes"] == ["Sox2", "Pax6", "Atoh1"]     # split on commas/whitespace, duplicates removed
        assert (save["user_id"], save["gctype"], save["label"], save["ldesc"]) == (1, "unweighted-list", "Hair cells",
                                                                                   "Markers")
        (record,) = result.logged("userhistory.add_record")
        assert record["user_id"] == 1 and record["gene_cart_share_id"] == save["share_id"]
        assert cart_files(www) == []

    def test_uploaded_unweighted(self, www):
        result = self.run(www, new_cart_upload_type="uploaded-unweighted",
                          files={"new_cart_file": ("genes.txt", b"Sox2\nPax6,Sox2\n")})
        assert result.json() == {"id": 7}
        assert self.saved(result)["genes"] == ["Sox2", "Pax6"]

    @pytest.mark.parametrize("upload_type", ["uploaded-unweighted", "uploaded-weighted"])
    def test_upload_without_a_file(self, www, upload_type):
        result = self.run(www, new_cart_upload_type=upload_type)
        self.assert_refused(result, www, "uploaded file")

    @pytest.mark.parametrize("filename, content", [
        ("markers.tab", WEIGHTED_TAB.encode()),
        ("markers.csv", WEIGHTED_TAB.replace("\t", ",").encode()),
        ("markers.xlsx", None),
    ])
    def test_uploaded_weighted(self, www, filename, content):
        if content is None:
            buffer = io.BytesIO()
            pd.read_csv(io.StringIO(WEIGHTED_TAB), sep="\t").to_excel(buffer, index=False)
            content = buffer.getvalue()
        result = self.run(www, new_cart_upload_type="uploaded-weighted", files={"new_cart_file": (filename, content)})
        assert_clean_exit(result)
        assert result.json() == {"id": 7}
        save = self.saved(result)
        assert save["gctype"] == "weighted-list" and save["user_id"] == 1
        share_id = save["share_id"]
        assert cart_files(www) == [f"cart.{share_id}.h5ad", f"cart.{share_id}.tab"]

        tab = pd.read_csv(www / "carts" / f"cart.{share_id}.tab", sep="\t")
        assert list(tab.columns) == ["id", "gene", "P1", "P2"]
        adata = anndata.read_h5ad(www / "carts" / f"cart.{share_id}.h5ad")
        assert list(adata.obs_names) == ["P1", "P2"]
        assert list(adata.var_names) == ["E1", "E2", "E3"] and list(adata.var["gene"]) == ["Sox2", "Pax6", "Atoh1"]
        assert adata.X.dtype == np.float64        # integer weights (P2) are stored as floats
        assert adata.X.tolist() == [[0.5, -1.5, 2.0], [1.0, 0.0, 1.0]]

    @pytest.mark.parametrize("filename, content, message", [
        ("markers.tab", "id\tgene\nE1\tSox2\n", "not valid"),                          # no weight columns
        ("markers.tab", "id\tgene\tP1\nE1\tSox2\t1\nE1\tPax6\t2\n", "not valid"),     # identifiers not unique
        ("markers.tab", "id\tgene\tP1\nE1\tSox2\thigh\nE2\tPax6\tlow\n", "not valid"),  # weights not numeric
        ("markers.tab", "id\tgene\tP1\n1\tSox2\t1\n2\tPax6\t2\n", "not valid"),       # identifiers not strings
        ("markers.txt", WEIGHTED_TAB, "Unsupported file type"),
    ])
    def test_invalid_weighted_upload(self, www, filename, content, message):
        result = self.run(www, new_cart_upload_type="uploaded-weighted",
                          files={"new_cart_file": (filename, content.encode())})
        self.assert_refused(result, www, message)

    @pytest.mark.parametrize("upload_type, message", [("labeled-list", "Not implemented"),
                                                      ("mystery", "Invalid upload type")])
    def test_unsupported_upload_type(self, www, upload_type, message):
        self.assert_refused(self.run(www, new_cart_upload_type=upload_type), www, message)

    def test_is_domain_is_not_read_from_the_form(self, www):
        # Only admins may make site-wide "domain" lists; this form has no such field and must not honor one
        result = self.run(www, new_cart_upload_type="pasted_genes", new_cart_pasted_genes="Sox2", is_domain="1")
        assert result.json() == {"id": 7}
        assert not self.saved(result)["is_domain"]


# ---------------------------------------------------------------- save_new_genecart_json.cgi

@pytest.mark.slow   # weighted lists are converted with scanpy in the CGI subprocess
class TestSaveNewGeneCartJson:
    DB = {"sessions": SESSIONS, "lastrowid": 9}
    UNWEIGHTED = {"label": "Hair cells", "organism_id": 1, "is_public": 0, "gctype": "unweighted-list",
                  "genes": [{"id": 11, "gene_symbol": "Sox2"}, {"id": 12, "gene_symbol": "Pax6"}]}
    WEIGHTED = {"label": "Patterns", "organism_id": 1, "is_public": 0, "gctype": "weighted-list",
                "weight_labels": ["P1", "P2"],
                "genes": [{"id": "E1", "gene_symbol": "Sox2", "weights": [0.5, 1]},
                          {"id": "E2", "gene_symbol": "Pax6", "weights": [-1.5, 0]}]}

    def run(self, www, payload):
        return run_json_cgi(copy_cgi(www, "save_new_genecart_json.cgi"), payload, db=self.DB)

    def test_unweighted(self, www):
        result = self.run(www, {"session_id": "owner", **self.UNWEIGHTED})
        assert_clean_exit(result)
        assert result.json() == {"id": 9}
        (save,) = result.logged("gene_cart.save")
        assert (save["user_id"], save["gctype"], save["label"], save["genes"]) == (1, "unweighted-list", "Hair cells",
                                                                                   ["Sox2", "Pax6"])
        assert cart_files(www) == []

    def test_weighted(self, www):
        result = self.run(www, {"session_id": "owner", **self.WEIGHTED})
        assert_clean_exit(result)
        assert result.json() == {"id": 9}
        (save,) = result.logged("gene_cart.save")
        share_id = save["share_id"]
        assert cart_files(www) == [f"cart.{share_id}.h5ad", f"cart.{share_id}.tab"]
        assert (www / "carts" / f"cart.{share_id}.tab").read_text() == (
            "identifier\tgene_symbol\tP1\tP2\nE1\tSox2\t0.5\t1\nE2\tPax6\t-1.5\t0\n")
        adata = anndata.read_h5ad(www / "carts" / f"cart.{share_id}.h5ad")
        assert list(adata.obs_names) == ["P1", "P2"] and list(adata.var_names) == ["E1", "E2"]

    @pytest.mark.parametrize("session", [None, "expired"])
    @pytest.mark.parametrize("cart", ["UNWEIGHTED", "WEIGHTED"])
    def test_anonymous_is_refused(self, www, session, cart):
        """Without a logged-in user nothing is saved: no gene list row with no owner, and no files in www/carts."""
        payload = dict(getattr(self, cart))
        if session:
            payload["session_id"] = session
        result = self.run(www, payload)
        assert result.status == 500
        assert "No logged-in user" in result.json()["message"]
        assert not result.logged("gene_cart.save")
        assert result.writes() == []
        assert cart_files(www) == []

    @pytest.mark.parametrize("session, requested, saved", [
        ("owner", 1, 0), ("owner", "1", 0), ("owner", True, 0),   # anyone else's request is ignored
        ("admin", 1, 1), ("admin", "1", 1), ("admin", True, 1),
        ("admin", 0, 0), ("admin", None, 0), ("admin", "x", 0),
    ])
    def test_only_admins_create_domain_lists(self, www, session, requested, saved):
        """Domain lists are shown to everyone on the site; update_from_json() copies is_domain from the request."""
        payload = {"session_id": session, **self.UNWEIGHTED}
        if requested is not None:
            payload["is_domain"] = requested
        result = run_json_cgi(copy_cgi(www, "save_new_genecart_json.cgi"), payload,
                              db={**self.DB, "sessions": {**SESSIONS, "admin": 3}, "users": {"3": {"is_admin": 1}}})
        assert_clean_exit(result)
        assert result.json() == {"id": 9}
        (save,) = result.logged("gene_cart.save")
        assert save["is_domain"] == saved

    def test_labeled_list_is_not_implemented(self, www):
        result = self.run(www, {"session_id": "owner", **self.UNWEIGHTED, "gctype": "labeled-list"})
        assert result.status == 500
        assert result.json() == {"error": "Not implemented"}
        assert not result.logged("gene_cart.save")


# ---------------------------------------------------------------- save_genecart_changes.cgi

class TestSaveGeneCartChanges:
    CART = {"id": 5, "share_id": "GC1", "user_id": 1, "gctype": "unweighted-list", "label": "Old",
            "organism_id": "1", "ldesc": "Old description", "is_public": 0}
    CHANGES = {"gc_id": "5", "visibility": "1", "title": "New", "organism_id": "2", "ldesc": "New description"}

    def run(self, cart=None, **form):
        db = {"sessions": SESSIONS, "gene_carts": [cart or self.CART]}
        return run_cgi("save_genecart_changes.cgi", db=db, form=form)

    def changes(self, result):
        return {e["attribute"]: e["value"] for e in result.logged("gene_cart.save_change")}

    def test_owner_saves_changes(self):
        result = self.run(session_id="owner", genes="Sox2, Pax6 Sox2", **self.CHANGES)
        body = result.json()
        assert body["success"] == 1
        assert json.loads(body["gene_cart"])["label"] == "New"     # the cart is returned as a JSON string
        assert self.changes(result) == {"is_public": 1, "label": "New", "organism_id": "2", "ldesc": "New description"}
        (members,) = result.logged("gene_cart.save_members")
        assert members["genes"] == ["Sox2", "Pax6"]

    def test_unchanged_fields_are_not_saved(self):
        result = self.run(session_id="owner", gc_id="5", visibility="0", title="Old", organism_id="1",
                          ldesc="Old description")
        assert result.json()["success"] == 1
        assert self.changes(result) == {}
        assert not result.logged("gene_cart.save_members")      # an empty gene list doesn't empty the cart

    def test_genes_are_ignored_for_weighted_lists(self):
        result = self.run(cart={**self.CART, "gctype": "weighted-list"}, session_id="owner", genes="Sox2",
                          **self.CHANGES)
        assert result.json()["success"] == 1
        assert not result.logged("gene_cart.save_members")

    @pytest.mark.parametrize("session, error", [(None, "No session id"), ("expired", "Invalid session id"),
                                                ("other", "You do not own this cart")])
    def test_non_owner_is_refused(self, session, error):
        form = {"genes": "Sox2", **self.CHANGES}
        if session:
            form["session_id"] = session
        result = self.run(**form)
        body = result.json()
        assert body["success"] == 0 and error in body["error"]
        assert not result.logged("gene_cart.save_change")
        assert not result.logged("gene_cart.save_members")
        assert result.writes() == []

    def test_unknown_cart(self):
        result = self.run(session_id="owner", **{**self.CHANGES, "gc_id": "99"})
        body = result.json()
        assert body["success"] == 0 and "Invalid gene cart id" in body["error"]
        assert not result.logged("gene_cart.save_change")


# ---------------------------------------------------------------- get_user_gene_carts.cgi

class TestGetUserGeneCarts:
    CARTS = [
        {"id": 1, "share_id": "DOM", "user_id": 9, "gctype": "unweighted-list", "is_domain": 1, "is_public": 1},
        {"id": 2, "share_id": "MINE", "user_id": 1, "gctype": "weighted-list", "recent": 1},
        {"id": 3, "share_id": "GROUP", "user_id": 3, "gctype": "unweighted-list", "group_user_ids": [1]},
        {"id": 4, "share_id": "PRIV", "user_id": 3, "gctype": "weighted-list"},
        {"id": 5, "share_id": "PUB", "user_id": 3, "gctype": "unweighted-list", "is_public": 1},
    ]
    KEYS = ["domain_carts", "group_carts", "public_carts", "recent_carts", "shared_carts", "user_carts"]

    def run(self, **query):
        return run_cgi("get_user_gene_carts.cgi", db={"sessions": SESSIONS, "gene_carts": self.CARTS}, query=query)

    def share_ids(self, body):
        return {key: [c["share_id"] for c in body[key]] for key in self.KEYS}

    def test_anonymous_sees_domain_public_and_shared(self):
        result = self.run(share_id="PRIV")
        assert self.share_ids(result.json()) == {
            "domain_carts": ["DOM"], "group_carts": [], "public_carts": ["DOM", "PUB"], "recent_carts": [],
            "shared_carts": ["PRIV"], "user_carts": []}
        assert not result.logged("gene_cart_collection.get_by_user")
        assert result.writes() == []

    def test_logged_in_user(self):
        body = self.run(session_id="owner").json()
        # A list can be in several categories (e.g. the user's recent lists are also in user_carts)
        assert self.share_ids(body) == {
            "domain_carts": ["DOM"], "group_carts": ["GROUP"], "public_carts": ["DOM", "PUB"],
            "recent_carts": ["MINE"], "shared_carts": [], "user_carts": ["MINE"]}

    @pytest.mark.parametrize("include_members, include_genes", [(None, True), ("1", True), ("0", False)])
    def test_include_members(self, include_members, include_genes):
        query = {"session_id": "owner"}
        if include_members is not None:
            query["include_members"] = include_members
        result = self.run(**query)
        assert result.json()["user_carts"]
        calls = [e for e in result.db_log if e["event"].startswith("gene_cart_collection.")]
        assert len(calls) == 6 and all(e["include_genes"] is include_genes for e in calls)

    def test_cart_type_filter(self):
        body = self.run(session_id="owner", share_id="PRIV", cart_type="weighted-list").json()
        assert self.share_ids(body) == {
            "domain_carts": [], "group_carts": [], "public_carts": [], "recent_carts": ["MINE"],
            "shared_carts": ["PRIV"], "user_carts": ["MINE"]}

    @pytest.mark.parametrize("cart_type", ["null", "mystery"])
    def test_unknown_cart_type_returns_everything(self, cart_type):
        result = self.run(session_id="owner", cart_type=cart_type)
        assert self.share_ids(result.json()) == self.share_ids(self.run(session_id="owner").json())

    def test_group_by_type(self):
        body = self.run(session_id="owner", group_by_type="true").json()
        assert list(body) == ["unweighted-list", "weighted-list", "labeled-list"]
        assert self.share_ids(body["weighted-list"]) == {
            "domain_carts": [], "group_carts": [], "public_carts": [], "recent_carts": ["MINE"],
            "shared_carts": [], "user_carts": ["MINE"]}
        assert self.share_ids(body["unweighted-list"])["group_carts"] == ["GROUP"]
        assert all(not carts for carts in body["labeled-list"].values())

    def test_group_by_type_false_is_not_grouped(self):
        body = self.run(session_id="owner", group_by_type="false").json()
        assert sorted(body) == self.KEYS


# ---------------------------------------------------------------- get_weighted_gene_cart_preview.cgi

@pytest.mark.slow   # weighted lists are converted with scanpy in the CGI subprocess
class TestGetWeightedGeneCartPreview:
    def run(self, www, share_id):
        query = {} if share_id is None else {"share_id": share_id}
        return run_cgi(copy_cgi(www, "get_weighted_gene_cart_preview.cgi"), query=query)

    def write_h5ad(self, www, share_id):
        adata = anndata.AnnData(X=np.array([[0.5, -1.5, 2.0], [1.0, 0.0, 1.0]]),
                                obs=pd.DataFrame(index=["P1", "P2"]),
                                var=pd.DataFrame({"gene": ["Sox2", "Pax6", "Atoh1"]}, index=["E1", "E2", "E3"]))
        adata.write_h5ad(www / "carts" / f"cart.{share_id}.h5ad")

    # Share IDs can be renamed by the owner to letters, digits, hyphens, underscores and periods (update_share_id.cgi)
    @pytest.mark.parametrize("share_id", ["ab12cd34", "my-markers_v1.2"])
    def test_preview(self, www, share_id):
        self.write_h5ad(www, share_id)
        result = self.run(www, share_id)
        assert_clean_exit(result)
        assert result.json() == {"preview_json": [], "success": 1, "num_genes": 3, "weights": ["P1", "P2"]}

    @pytest.mark.parametrize("share_id", ["NOPE", "", None, "../ab12cd34"])
    def test_missing_cart(self, www, share_id):
        self.write_h5ad(www, "ab12cd34")
        result = self.run(www, share_id)
        assert_clean_exit(result)
        body = result.json()
        assert body["success"] == -1 and body["num_genes"] == -1
        assert body["message"].startswith("Unable to open file")


# ---------------------------------------------------------------- get_unweighted_gene_cart_preview.cgi

class TestGetUnweightedGeneCartPreview:
    CART = {"id": 5, "share_id": "ab12cd34", "user_id": 1, "gctype": "unweighted-list", "organism_id": 1,
            "genes": ["Sox2", "Pax6"]}
    GENES = [{"gene_symbol": "Sox2", "ensembl_id": "ENSMUSG01", "product": "SRY-box 2", "organism_id": 1},
             {"gene_symbol": "Pax6", "ensembl_id": "ENSMUSG02", "product": "paired box 6", "organism_id": 1}]

    def run(self, share_id, cart=None):
        query = {} if share_id is None else {"share_id": share_id}
        return run_cgi("get_unweighted_gene_cart_preview.cgi", query=query,
                       db={"gene_carts": [cart or self.CART], "genes": self.GENES})

    @pytest.mark.parametrize("share_id", ["ab12cd34", "my-markers_v1.2"])
    def test_preview(self, share_id):
        result = self.run(share_id, cart={**self.CART, "share_id": share_id})
        assert result.json() == {"success": 1, "gene_info": {
            "ENSMUSG01": {"gene_symbol": "Sox2", "product": "SRY-box 2"},
            "ENSMUSG02": {"gene_symbol": "Pax6", "product": "paired box 6"}}}
        (search,) = result.logged("gene_collection.get_by_gene_symbol")
        assert search["gene_symbol"] == "Sox2 Pax6" and search["exact"] is True

    @pytest.mark.parametrize("share_id, message", [("NOPE", "Invalid share ID"), ("", "No share ID provided"),
                                                   (None, "No share ID provided"), ("/", "No share ID provided")])
    def test_missing_cart(self, share_id, message):
        result = self.run(share_id)
        assert_clean_exit(result)
        assert result.json() == {"success": 0, "message": message}
        assert not result.logged("gene_collection.get_by_gene_symbol")


# ---------------------------------------------------------------- get_pattern_element_list.cgi

class TestGetPatternElementList:
    def run(self, **query):
        return run_cgi("get_pattern_element_list.cgi", query=query)

    def test_unweighted_list(self):
        body = self.run(source_id="anything", scope="unweighted-list").json()
        assert body == [{"label": "unweighted", "top_up": "", "top_down": "", "binary": True}]

    def test_weighted_list(self, tmp_cart):
        share_id = tmp_cart(WEIGHTED_TAB)
        result = self.run(source_id=share_id, scope="weighted-list")
        assert_clean_exit(result)
        assert result.json() == [
            {"label": "P1", "top_up": "Atoh1,Sox2,Pax6", "top_down": "Pax6,Sox2,Atoh1", "binary": False},
            {"label": "P2", "top_up": "Sox2,Atoh1,Pax6", "top_down": "", "binary": True},   # no negative weights
        ]

    @pytest.mark.parametrize("query", [{"source_id": "pytest-no-such-cart"}, {}, {"scope": "weighted-list"}])
    def test_missing_cart_is_empty(self, query):
        result = self.run(**query)
        assert_clean_exit(result)
        assert result.json() == []
