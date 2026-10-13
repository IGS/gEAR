"""
Single-cell workbench step CGIs (h5ad_*), PCs lookup, stored analyses and upload-time obs dtype review.

Step CGIs run against fakes/cgi_wb_analysis.py: analyses live under $GEAR_TEST_ANALYSIS_DIR and a
primary analysis is a small synthetic AnnData (two louvain groups, three mt- genes). Stored-analysis
CGIs run against fakes/cgi_wb_real_analysis.py, the real gear.analysis moved into a temporary tree.
"""

import importlib.machinery
import importlib.util
import json
import sys
import uuid
from pathlib import Path

import pytest
from helpers.cgi_harness import CGI_DIR, run_cgi

# Every workbench step CGI imports scanpy in its subprocess
pytestmark = pytest.mark.slow

FAKES = Path(__file__).resolve().parent / "fakes"
STEP_ANALYSIS = {"gear.analysis": FAKES / "cgi_wb_analysis.py"}
REAL_ANALYSIS = {"gear.analysis": FAKES / "cgi_wb_real_analysis.py"}
DB = {"sessions": {"owner": 1, "other": 2}, "datasets": [{"id": "DS1", "owner_id": 1}]}


def files_under(path):
    return sorted(str(p.relative_to(path)) for p in Path(path).rglob("*") if p.is_file())


@pytest.fixture
def adata_builder():
    pytest.importorskip("anndata")
    from fakes.cgi_wb_analysis import build_adata
    return build_adata


def write_analysis(tmp_path, adata, *parts):
    """Write adata as <tmp>/<parts...>/DS1.h5ad (the fake's layout) and return the path."""
    path = tmp_path.joinpath(*parts, "DS1.h5ad")
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    return path


def run_step(tmp_path, script, db=DB, **query):
    return run_cgi(script, db=db, query=query, fake_modules=STEP_ANALYSIS,
                   env={"GEAR_TEST_ANALYSIS_DIR": str(tmp_path)})


def unsaved_dir(tmp_path, analysis_id="A1"):
    return tmp_path / "user_unsaved" / "owner" / analysis_id


# Minimal valid parameters for each step, from the primary analysis
STEP_QUERIES = {
    "h5ad_generate_pca.cgi": {"compute_pca": "true"},
    "h5ad_identify_variable_genes.cgi": {"norm_counts_per_cell": "1e4", "flavor": "seurat", "n_top_genes": "5",
                                         "min_mean": "0.0125", "max_mean": "3", "min_dispersion": "0.5",
                                         "regress_out": "false", "scale_unit_variance": "false",
                                         "save_dataset": "1"},
    "h5ad_qc_by_mito.cgi": {"genes_prefix": "mt-", "filter_mito_perc": "", "filter_mito_count": "",
                            "save_dataset": "1"},
    "h5ad_preview_primary_filter.cgi": {},
    "h5ad_compare_genes.cgi": {"query_cluster": "g1", "reference_cluster": "g2", "n_genes": "3",
                               "method": "t-test", "corr_method": "benjamini-hochberg"},
}


# ---------------------------------------------------------------- errors shared by every step

class TestStepRefusals:
    @pytest.mark.parametrize("script", STEP_QUERIES)
    def test_unknown_dataset(self, tmp_path, script):
        query = {"dataset_id": "NOPE", "analysis_id": "A1", "analysis_type": "primary", "session_id": "owner",
                 **STEP_QUERIES[script]}
        result = run_step(tmp_path, script, **query)
        assert result.json()["success"] == 0, result.stderr[-2000:]
        assert files_under(tmp_path) == []
        assert result.writes() == []

    @pytest.mark.parametrize("script", STEP_QUERIES)
    def test_other_users_saved_analysis_is_unavailable(self, tmp_path, adata_builder, script):
        # user_saved analyses are found under the session's user, so another user can't reach the owner's
        owned = write_analysis(tmp_path, adata_builder(), "user_saved", "1", "A1")
        before = owned.read_bytes()
        query = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "user_saved", "session_id": "other",
                 **STEP_QUERIES[script]}
        result = run_step(tmp_path, script, **query)
        assert result.json()["success"] == 0, result.stderr[-2000:]
        assert files_under(tmp_path) == ["user_saved/1/A1/DS1.h5ad"]
        assert owned.read_bytes() == before
        assert result.writes() == []


# ---------------------------------------------------------------- h5ad_preview_primary_filter.cgi

class TestPreviewPrimaryFilter:
    def run(self, tmp_path, db=DB):
        return run_step(tmp_path, "h5ad_preview_primary_filter.cgi", db=db, dataset_id="DS1",
                        analysis_id="DS1", analysis_type="primary", session_id="owner")

    def test_success_when_prelim_violin_exists(self, tmp_path):
        (tmp_path / "primary").mkdir()
        (tmp_path / "primary" / "DS1.prelim_violin.png").write_bytes(b"png")
        assert self.run(tmp_path).json() == {"success": 1, "is_spatial": 0}

    def test_no_violin_is_not_success(self, tmp_path):
        assert self.run(tmp_path).json() == {"success": 0, "is_spatial": 0}

    def test_reports_spatial_dataset(self, tmp_path):
        db = {**DB, "datasets": [{"id": "DS1", "owner_id": 1, "dtype": "spatial"}]}
        assert self.run(tmp_path, db=db).json() == {"success": 0, "is_spatial": 1}


# ---------------------------------------------------------------- h5ad_qc_by_mito.cgi

class TestQcByMito:
    @pytest.fixture(autouse=True)
    def needs_scanpy(self):
        pytest.importorskip("scanpy")

    def run(self, tmp_path, **query):
        query = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "primary", "session_id": "owner",
                 **STEP_QUERIES["h5ad_qc_by_mito.cgi"], **query}
        result = run_step(tmp_path, "h5ad_qc_by_mito.cgi", **query)
        assert result.returncode == 0, result.stderr[-3000:]
        return result

    @staticmethod
    def expected(adata):
        import numpy as np
        mito = adata.var.gene_symbol.str.lower().str.startswith("mt-").to_numpy()
        return np.asarray(adata.X[:, mito].sum(axis=1) / adata.X.sum(axis=1)).ravel()

    def test_filters_and_saves_to_unsaved_analysis(self, tmp_path, adata_builder):
        import anndata
        import numpy as np
        adata = adata_builder()
        percent = self.expected(adata)
        cutoff = float(np.median(percent))
        max_genes = int(adata.obs["n_genes"].max())
        keep = (adata.obs["n_genes"].to_numpy() < max_genes) & (percent < cutoff)

        # The prefix is matched case-insensitively
        body = self.run(tmp_path, genes_prefix="MT-", filter_mito_perc=repr(cutoff),
                        filter_mito_count=str(max_genes)).json()
        assert body == {"success": 1, "n_obs": int(keep.sum()), "n_genes": adata.n_vars}

        # The primary analysis is never written; the result goes to the user's unsaved analysis
        saved = anndata.read_h5ad(unsaved_dir(tmp_path) / "DS1.h5ad")
        assert np.allclose(saved.obs["percent_mito"], percent[keep])
        assert not (tmp_path / "primary").exists()
        assert sorted(p.name for p in (unsaved_dir(tmp_path) / "figures").iterdir()) == [
            "scatter_n_genes.png", "scatter_percent_mito.png", "violin_qc_by_mito.png"]

    def test_preview_does_not_save(self, tmp_path):
        assert self.run(tmp_path, save_dataset="0").json() == {"success": 1}
        assert not (unsaved_dir(tmp_path) / "DS1.h5ad").exists()


# ---------------------------------------------------------------- h5ad_identify_variable_genes.cgi

class TestIdentifyVariableGenes:
    @pytest.fixture(autouse=True)
    def needs_scanpy(self):
        pytest.importorskip("scanpy")

    def run(self, tmp_path, **query):
        query = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "primary", "session_id": "owner",
                 **STEP_QUERIES["h5ad_identify_variable_genes.cgi"], **query}
        result = run_step(tmp_path, "h5ad_identify_variable_genes.cgi", **query)
        assert result.returncode == 0, result.stderr[-3000:]
        return result.json()

    def test_top_n_genes_saved(self, tmp_path):
        import anndata
        body = self.run(tmp_path)
        assert body["success"] == 1 and (body["n_obs"], body["n_genes"]) == (60, 5)
        top = body["top_genes"].split(", ")
        assert len(top) == 5
        saved = anndata.read_h5ad(unsaved_dir(tmp_path) / "DS1.h5ad")
        # Every gene is kept in the file, with the highly variable ones flagged
        assert saved.n_vars == 30 and int(saved.var["highly_variable"].sum()) == 5
        assert set(top) <= set(saved.var.gene_symbol[saved.var["highly_variable"]])
        assert saved.raw is not None
        assert (unsaved_dir(tmp_path) / "figures" / "filter_genes_dispersion.png").exists()

    def test_regress_out_and_scale(self, tmp_path, adata_builder):
        import anndata
        import numpy as np
        # Continues an unsaved analysis whose QC step stored n_counts and percent_mito
        adata = adata_builder()
        adata.obs["n_counts"] = adata.X.sum(axis=1)
        adata.obs["percent_mito"] = adata.X[:, -3:].sum(axis=1) / adata.obs["n_counts"]
        write_analysis(tmp_path, adata, "user_unsaved", "owner", "A1")
        body = self.run(tmp_path, analysis_type="user_unsaved", regress_out="true", scale_unit_variance="true")
        assert body["success"] == 1
        saved = anndata.read_h5ad(unsaved_dir(tmp_path) / "DS1.h5ad")
        assert np.allclose(saved.X.mean(axis=0), 0, atol=1e-4) and saved.X.max() <= 10

    def test_dispersion_cutoffs_without_saving(self, tmp_path):
        body = self.run(tmp_path, n_top_genes="", save_dataset="0")
        assert body["success"] == 1 and body["n_obs"] == 60
        assert 0 <= body["n_genes"] <= 30
        assert not (unsaved_dir(tmp_path) / "DS1.h5ad").exists()


# ---------------------------------------------------------------- h5ad_generate_pca.cgi

class TestGeneratePca:
    @pytest.fixture(autouse=True)
    def needs_scanpy(self):
        pytest.importorskip("scanpy")

    def run(self, tmp_path, **query):
        query = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "primary", "session_id": "owner",
                 **query}
        result = run_step(tmp_path, "h5ad_generate_pca.cgi", **query)
        assert result.returncode == 0, result.stderr[-3000:]
        return result.json()

    def test_compute_then_replot_from_saved_pca(self, tmp_path):
        import anndata
        # Gene names are matched case-insensitively; unknown ones are ignored
        assert self.run(tmp_path, compute_pca="true", genes_to_color="g0, G1,NOPE") == \
            {"success": 1, "missing_gene": ""}
        saved = anndata.read_h5ad(unsaved_dir(tmp_path) / "DS1.h5ad")
        assert "X_pca" in saved.obsm and "PCs" in saved.varm
        figures = unsaved_dir(tmp_path) / "figures"
        assert sorted(p.name for p in figures.iterdir()) == ["pca.png", "pca_variance_ratio.png"]

        # Without compute_pca the stored PCA is plotted again
        (figures / "pca.png").unlink()
        assert self.run(tmp_path, analysis_type="user_unsaved") == {"success": 1, "missing_gene": ""}
        assert (figures / "pca.png").exists()


def _load_cgi(name):
    """Import a CGI as a module, restoring sys.stdout (step CGIs redirect it to /dev/null on import)."""
    stdout = sys.stdout
    try:
        loader = importlib.machinery.SourceFileLoader(name.replace(".", "_"), str(CGI_DIR / name))
        module = importlib.util.module_from_spec(importlib.util.spec_from_loader(loader.name, loader))
        loader.exec_module(module)
    finally:
        if sys.stdout is not stdout:
            sys.stdout.close()
            sys.stdout = stdout
    return module


class TestPcaHelpers:
    @pytest.fixture(scope="class")
    def pca(self):
        pytest.importorskip("scanpy")
        return _load_cgi("h5ad_generate_pca.cgi")

    def test_normalize_genes_to_color(self, pca):
        genes = ["Sox2", "Pax6", "mt-Co1"]
        assert pca.normalize_genes_to_color(genes, ["pax6", "SOX2", "Missing"]) == ["Pax6", "Sox2"]
        assert pca.normalize_genes_to_color(genes, []) == []

    @pytest.mark.parametrize("num_pcs", [1, 2, 5])
    def test_broken_stick_model(self, pca, num_pcs):
        # Expected share of component k (1-based) is (1/p) * sum_{i=k..p} 1/i, scaled to the PCs' total
        expected = [sum(1 / i for i in range(k, num_pcs + 1)) / num_pcs for k in range(1, num_pcs + 1)]
        ratios = pca.create_broken_stick_model(num_pcs, 0.5)
        assert ratios == pytest.approx([0.5 * e for e in expected])
        assert sum(ratios) == pytest.approx(0.5)
        assert ratios == sorted(ratios, reverse=True)


# ---------------------------------------------------------------- get_PCs_from_anndata.cgi

class TestGetPCs:
    def run(self, tmp_path, **query):
        query = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "user_unsaved", "session_id": "owner",
                 **query}
        return run_step(tmp_path, "get_PCs_from_anndata.cgi", **query)

    def test_returns_loadings_by_gene(self, tmp_path, adata_builder):
        import numpy as np
        adata = adata_builder()
        adata.varm["PCs"] = np.arange(adata.n_vars * 2, dtype=float).reshape(adata.n_vars, 2)
        write_analysis(tmp_path, adata, "user_unsaved", "owner", "A1")
        body = self.run(tmp_path).json()
        assert body["success"] == 1
        assert body["pc_data"]["columns"] == ["PC1", "PC2"]
        assert body["pc_data"]["index"] == adata.var.index.tolist()
        assert body["pc_data"]["data"][1] == [2.0, 3.0]
        assert body["gene_symbols"] == adata.var.gene_symbol.tolist()

    @pytest.mark.parametrize("query, message", [
        ({"dataset_id": ""}, "No dataset_id provided"),
        ({"dataset_id": "NOPE"}, "No dataset found with that ID"),
        ({}, "Analysis for this dataset is unavailable."),                      # no stored analysis
        ({"analysis_id": "", "analysis_type": "primary"}, "Analysis for this dataset is unavailable."),
        ({"analysis_type": "primary"}, "PCs not found in AnnData object"),     # PCA never run
    ])
    def test_errors(self, tmp_path, query, message):
        pytest.importorskip("anndata")
        assert self.run(tmp_path, **query).json() == {"success": 0, "message": message, "pc_data": None,
                                                     "gene_symbols": None}


# ---------------------------------------------------------------- h5ad_compare_genes.cgi

class TestCompareGenes:
    @pytest.fixture(autouse=True)
    def needs_scanpy(self):
        pytest.importorskip("scanpy")

    def run(self, tmp_path, **query):
        query = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "primary", "session_id": "owner",
                 **STEP_QUERIES["h5ad_compare_genes.cgi"], **query}
        result = run_step(tmp_path, "h5ad_compare_genes.cgi", **query)
        assert result.returncode == 0, result.stderr[-3000:]
        return result.json()

    @staticmethod
    def table(body, key):
        return json.loads(body[key])

    def test_query_against_reference_cluster(self, tmp_path):
        body = self.run(tmp_path)
        assert body["success"] == 1 and body["cluster_label"] == "louvain"
        forward, reverse = self.table(body, "table_json_f"), self.table(body, "table_json_r")
        assert forward["columns"] == ["names", "log2FC", "P-value", "FDR"]
        # n_genes rows, led by the genes made higher in g1
        assert len(forward["data"]) == 3 and len(reverse["data"]) == 3
        assert {row[0] for row in forward["data"]} <= {"G0", "G1", "G2", "G3", "G4"}
        assert {row[0] for row in reverse["data"]}.isdisjoint({"G0", "G1", "G2", "G3", "G4"})
        figures = sorted(p.name for p in (unsaved_dir(tmp_path) / "figures").iterdir())
        assert figures == ["rank_genes_groups_louvain_comp_ranked.png", "rank_genes_groups_louvain_comp_ranked_rev.png",
                           "rank_genes_groups_louvain_g1_comp_violin.png",
                           "rank_genes_groups_louvain_g2_comp_violin_rev.png"]

    def test_renamed_groups_against_all_others(self, tmp_path):
        body = self.run(tmp_path, group_labels='["first", "second"]', query_cluster="first",
                        reference_cluster="all-reference-clusters", n_genes="2")
        assert body["success"] == 1
        assert len(self.table(body, "table_json_f")["data"]) == 2
        figures = sorted(p.name for p in (unsaved_dir(tmp_path) / "figures").iterdir())
        assert figures == ["rank_genes_groups_louvain_comp_ranked.png",
                           "rank_genes_groups_louvain_first_comp_violin.png"]

    def test_logistic_regression_has_no_tables(self, tmp_path):
        body = self.run(tmp_path, method="logreg")
        assert body["success"] == 1
        assert body["table_json_f"] is None and body["table_json_r"] is None

    def test_logistic_regression_needs_a_reference_cluster(self, tmp_path):
        body = self.run(tmp_path, method="logreg", reference_cluster="all-reference-clusters")
        assert body == {"success": 0, "cluster_label": "", "table_json_f": None, "table_json_r": None,
                        "error": "Logistic regression method requires a specific reference cluster."}
        assert files_under(tmp_path) == []


# ---------------------------------------------------------------- stored analyses (real gear.analysis)

U1, U2, U3, P1, P2 = (str(uuid.UUID(int=i)) for i in range(1, 6))

# The fake geardb's users have no is_curator flag; this module (written next to the analysis tree and
#  loaded with fake_modules) reads one from the spec: "users": {"3": {"is_curator": 1}}
CURATOR_PATCH = """
import sys
_geardb = sys.modules["geardb"]
_geardb._User.is_curator = property(
    lambda self: _geardb._SPEC.get("users", {}).get(str(self.id), {}).get("is_curator", 0))
"""


class StoredAnalyses:
    """Pipeline files in a temporary tree, read and written through the real gear.analysis."""

    @pytest.fixture
    def root(self, tmp_path):
        primary = tmp_path / "www" / "datasets"
        primary.mkdir(parents=True)
        (primary / "DS1.h5ad").write_bytes(b"")
        return tmp_path

    @staticmethod
    def db(root):
        return {**DB, "sessions": {**DB["sessions"], "curator": 3}, "users": {"3": {"is_curator": 1}},
                "datasets": [{"id": "DS1", "owner_id": 1, "h5ad_path": str(root / "www/datasets/DS1.h5ad")}]}

    @staticmethod
    def area(root, analysis_type, analysis_id, owner=None):
        return {
            "primary": root / "www" / "datasets",
            "public": root / "www" / "analyses" / "by_dataset" / "DS1" / str(analysis_id),
            "user_saved": root / "www" / "analyses" / "by_user" / str(owner) / "DS1" / str(analysis_id),
            "user_unsaved": root / "tmp" / str(owner) / "DS1" / str(analysis_id),
        }[analysis_type]

    def store(self, root, analysis_type, analysis_id, owner=None, h5ad=True, **extra):
        """Write a pipeline JSON (and an empty h5ad) for an analysis; returns the analysis directory."""
        path = self.area(root, analysis_type, analysis_id, owner)
        path.mkdir(parents=True, exist_ok=True)
        session = owner if analysis_type == "user_unsaved" else "owner"
        pipeline = {"id": analysis_id, "label": f"{analysis_type} {analysis_id}", "type": analysis_type,
                    "dataset_id": "DS1", "user_session_id": session, **extra}
        (path / "DS1.pipeline.json").write_text(json.dumps(pipeline))
        if h5ad:
            (path / "DS1.h5ad").write_bytes(b"")
        return path

    def run(self, root, script, **query):
        patch = root / "fake_curators.py"
        patch.write_text(CURATOR_PATCH)
        return run_cgi(script, db=self.db(root), query=query, fake_modules={"fake_curators": patch, **REAL_ANALYSIS},
                       env={"GEAR_TEST_ANALYSIS_ROOT": str(root)})


class TestGetStoredAnalysis(StoredAnalyses):
    def get(self, root, **query):
        return self.run(root, "get_stored_analysis.cgi", **{"dataset_id": "DS1", "session_id": "owner", **query})

    def test_returns_pipeline_json(self, root):
        self.store(root, "user_unsaved", U1, owner="owner", group_labels={"0": "a", "1": "b"})
        body = self.get(root, analysis_id=U1, analysis_type="user_unsaved").json()
        assert body == {"id": U1, "label": f"user_unsaved {U1}", "type": "user_unsaved", "dataset_id": "DS1",
                        "group_labels": ["a", "b"], "success": 1, "vetting": "owner"}

    @pytest.mark.parametrize("publisher, viewer, expected", [
        ("owner", "owner", "owner"),
        ("owner", "other", "community"),
        ("owner", None, "community"),
        ("curator", "other", "gear"),
    ])
    def test_public_analysis_vetting_comes_from_its_publisher(self, root, publisher, viewer, expected):
        """Vetting is derived from who published the analysis; a vetting saved in the pipeline file is ignored."""
        self.store(root, "public", P1, analysis_session_id=publisher, vetting="gear" if expected != "gear" else "owner")
        query = {"analysis_id": P1, "analysis_type": "public"}
        query["session_id"] = viewer
        if viewer is None:
            del query["session_id"]
        body = self.run(root, "get_stored_analysis.cgi", dataset_id="DS1", **query).json()
        assert body["success"] == 1
        assert body["vetting"] == expected

    def test_saved_analysis_vetting_ignores_value_in_file(self, root):
        self.store(root, "user_saved", U1, owner=1, vetting="gear")
        body = self.get(root, analysis_id=U1, analysis_type="user_saved").json()
        assert body["vetting"] == "owner"

    def test_public_analysis_does_not_reveal_the_publishers_session(self, root):
        # Anyone can load a public analysis; its pipeline file holds the publishing user's login session
        self.store(root, "public", P1, analysis_session_id="owner")
        body = self.get(root, session_id="other", analysis_id=P1, analysis_type="public").json()
        assert body["success"] == 1
        assert "analysis_session_id" not in body and "user_session_id" not in body

    def test_primary_analysis(self, root):
        self.store(root, "primary", "DS1", h5ad=False)
        body = self.get(root).json()
        assert body["success"] == 1 and body["type"] == "primary" and body["vetting"] == "gear"

    def test_missing_pipeline_file_is_one_json_error(self, root):
        # This used to print a second Content-Type header inside the body
        path = self.store(root, "user_unsaved", U1, owner="owner")
        (path / "DS1.pipeline.json").unlink()
        result = self.get(root, analysis_id=U1, analysis_type="user_unsaved")
        assert result.json() == {"success": 0, "error": "Analysis config file not found"}

    def test_unknown_dataset(self, root):
        assert self.get(root, dataset_id="NOPE").json() == {"success": 0, "error": "Dataset not found"}

    def test_other_users_saved_analysis_is_not_found(self, root):
        self.store(root, "user_saved", U1, owner=1)
        result = self.get(root, session_id="other", analysis_id=U1, analysis_type="user_saved")
        assert result.json() == {"success": 0, "error": "Analysis could not be loaded"}

    def test_path_in_analysis_id_cannot_reach_other_users_analysis(self, root):
        self.store(root, "user_saved", U1, owner=1)
        self.store(root, "user_saved", U2, owner=2)     # so the requesting user's own directory exists
        result = self.get(root, session_id="other", analysis_id=f"../../1/DS1/{U1}", analysis_type="user_saved")
        assert result.json()["success"] == 0


class TestGetStoredAnalysisList(StoredAnalyses):
    @pytest.fixture
    def tree(self, root):
        self.store(root, "primary", "DS1", h5ad=False)
        self.store(root, "public", P1)
        self.store(root, "user_saved", U1, owner=1)
        self.store(root, "user_saved", U2, owner=2)
        self.store(root, "user_unsaved", U3, owner="owner")
        return root

    def listing(self, root, **query):
        body = self.run(root, "get_stored_analysis_list.cgi", **{"dataset_id": "DS1", **query}).json()
        return {atype: [(a["id"], a.get("vetting")) for a in analyses] for atype, analyses in body.items()}

    def test_owner_sees_own_analyses(self, tree):
        # The owner comes from the session that wrote each pipeline file, or from the area it is in
        assert self.listing(tree, session_id="owner") == {
            "primary": [("DS1", "gear")], "public": [(P1, "owner")],
            "user_saved": [(U1, "owner")], "user_unsaved": [(U3, "owner")]}

    def test_other_user_sees_only_their_own(self, tree):
        # U2's pipeline was written under user 1's session, but it is in user 2's saved area
        assert self.listing(tree, session_id="other") == {
            "primary": [("DS1", "gear")], "public": [(P1, "community")], "user_saved": [(U2, "owner")],
            "user_unsaved": []}

    @pytest.mark.parametrize("session_id", ["owner", "other", "curator"])
    def test_public_analysis_by_a_curator_is_gear_vetted(self, tree, session_id):
        self.store(tree, "public", P2, user_session_id="curator")
        public = dict(self.listing(tree, session_id=session_id)["public"])
        assert public[P2] == "gear"

    def test_curator_viewing_another_users_analysis_does_not_make_it_gear_vetted(self, tree):
        # discover_vetting() used to check whether the viewer, not the owner, was a curator
        assert self.listing(tree, session_id="curator")["public"] == [(P1, "community")]

    def test_anonymous_sees_primary_and_public(self, tree):
        assert self.listing(tree) == {"primary": [("DS1", None)], "public": [(P1, None)],
                                      "user_saved": [], "user_unsaved": []}

    @pytest.mark.parametrize("session_id", [None, "other"])
    def test_publishers_session_is_not_listed(self, tree, session_id):
        query = {"session_id": session_id} if session_id else {}
        body = self.run(tree, "get_stored_analysis_list.cgi", dataset_id="DS1", **query).json()
        assert all("session_id" not in ana for analyses in body.values() for ana in analyses)

    def test_unknown_dataset(self, root):
        result = self.run(root, "get_stored_analysis_list.cgi", dataset_id="NOPE", session_id="owner")
        assert result.json() == {"success": 0, "error": "No dataset found with that ID"}


class TestSaveDatasetAnalysis(StoredAnalyses):
    STATE = '{"id": "x", "label": "Mine"}'

    def save(self, root, **query):
        query = {"dataset_id": "DS1", "session_id": "owner", "analysis_id": U1, "analysis_type": "user_unsaved",
                 "analysis_vetting": "undefined", "label": "Mine", "state": self.STATE, **query}
        return self.run(root, "save_dataset_analysis.cgi", **query)

    @pytest.mark.parametrize("analysis_type, owner", [("user_unsaved", "owner"), ("user_saved", 1)])
    def test_saves_pipeline(self, root, analysis_type, owner):
        path = self.area(root, analysis_type, U1, owner)
        path.mkdir(parents=True)
        assert self.save(root, analysis_type=analysis_type).json() == {"success": 1, "error": ""}
        assert (path / "DS1.pipeline.json").read_text() == self.STATE

    def test_write_failure_is_reported(self, root):
        body = self.save(root).json()        # the analysis directory doesn't exist
        assert body == {"success": 0,
                        "error": "An error occurred when attempting to write analysis pipeline path file"}

    def test_not_logged_in(self, root):
        self.area(root, "user_unsaved", U1, "nobody").mkdir(parents=True)
        assert self.save(root, session_id="nobody").json() == {"success": 0, "error": "User not found"}
        assert files_under(root / "tmp") == []

    def test_unknown_dataset(self, root):
        assert self.save(root, dataset_id="NOPE").json()["success"] == 0

    @pytest.mark.parametrize("analysis_type, analysis_id", [
        ("primary", "DS1"),                         # the dataset's own pipeline, shared with everyone
        ("public", P1),
        ("user_unsaved", "../../../www/datasets"),  # resolves to the primary pipeline
        ("user_unsaved", "not-a-uuid"),
        ("bogus", U1),
    ])
    def test_shared_or_odd_targets_are_refused(self, root, analysis_type, analysis_id):
        primary = self.store(root, "primary", "DS1", h5ad=False) / "DS1.pipeline.json"
        public = self.store(root, "public", P1) / "DS1.pipeline.json"
        self.area(root, "user_unsaved", "not-a-uuid", "owner").mkdir(parents=True)
        before = {primary: primary.read_text(), public: public.read_text()}
        body = self.save(root, analysis_type=analysis_type, analysis_id=analysis_id).json()
        assert body == {"success": 0, "error": "Only your own saved or unsaved analyses can be saved"}
        assert {path: path.read_text() for path in before} == before
        assert files_under(root / "tmp") == []


# ---------------------------------------------------------------- upload-time obs dtype review

QUESTIONABLE = {"replicate": {"current": "continuous", "suggested": "categorical"}}


@pytest.fixture
def staged(upload_session):
    """A staged upload <session>/<share>/ with metadata.json and an h5ad whose 'replicate' column is numeric."""
    session_id, path = upload_session
    share = path / "pytestshare"
    share.mkdir()
    (share / "metadata.json").write_text(json.dumps({"dataset_format": "single-cell-rnaseq",
                                                     "questionable_obs_columns": QUESTIONABLE}))
    return session_id, share


def stage_h5ad(share):
    anndata = pytest.importorskip("anndata")
    import numpy as np
    import pandas as pd
    obs = pd.DataFrame({"replicate": [1, 2, 1, 2], "dose": ["0.5", "1", "2", "x"]},
                       index=[f"c{i}" for i in range(4)])
    anndata.AnnData(X=np.ones((4, 3), dtype=np.float32), obs=obs).write_h5ad(share / f"{share.name}.h5ad")


class TestGetAmbiguousObsColumns:
    def run(self, session_id, share_uid="pytestshare", logged_in=True):
        db = {"sessions": {session_id: 1}} if logged_in else {}
        return run_cgi("get_ambiguous_obs_columns.cgi", db=db,
                       form={"session_id": session_id, "share_uid": share_uid})

    def test_returns_flagged_columns(self, staged):
        session_id, _ = staged
        assert self.run(session_id).json() == {"success": 1, "message": "", "questionable_columns": QUESTIONABLE,
                                                "reviewed": False}

    @pytest.mark.parametrize("share_uid, logged_in, message", [
        ("", True, "session_id and share_uid are required."),
        ("pytestshare", False, "Invalid session_id."),
        ("other", True, "No metadata JSON file found for this dataset."),
        ("../x/pytestshare", True, "No metadata JSON file found for this dataset."),  # sanitized to x_pytestshare
    ])
    def test_refusals(self, staged, share_uid, logged_in, message):
        session_id, _ = staged
        body = self.run(session_id, share_uid, logged_in).json()
        assert body["success"] == 0 and body["message"] == message

    def test_unreadable_metadata(self, staged):
        session_id, share = staged
        (share / "metadata.json").write_text("{not json")
        body = self.run(session_id).json()
        assert body["success"] == 0 and body["message"].startswith("Error reading metadata:")


class TestApplyObsDtypeChoices:
    def run(self, session_id, choices, share_uid="pytestshare", logged_in=True):
        db = {"sessions": {session_id: 1}} if logged_in else {}
        form = {"session_id": session_id, "share_uid": share_uid}
        if choices is not None:
            form["choices"] = choices if isinstance(choices, str) else json.dumps(choices)
        return run_cgi("apply_obs_dtype_choices.cgi", db=db, form=form)

    def test_applies_choices_and_marks_reviewed(self, staged):
        session_id, share = staged
        stage_h5ad(share)
        body = self.run(session_id, {"replicate": "categorical", "dose": "continuous"}).json()
        assert body == {"success": 1, "message": "Column types updated."}

        import anndata
        import pandas as pd
        obs = anndata.read_h5ad(share / "pytestshare.h5ad").obs
        assert isinstance(obs["replicate"].dtype, pd.CategoricalDtype)
        assert list(obs["replicate"]) == ["1", "2", "1", "2"]
        assert obs["dose"].iloc[:3].tolist() == [0.5, 1.0, 2.0] and pd.isna(obs["dose"].iloc[3])
        metadata = json.loads((share / "metadata.json").read_text())
        assert metadata["obs_dtype_reviewed"] is True and metadata["questionable_obs_columns"] == {}

    @pytest.mark.parametrize("choices, logged_in, message", [
        (None, True, "session_id, share_uid, and choices are required."),
        ({"replicate": "categorical"}, False, "Invalid session_id."),
        ("{bad", True, "Invalid choices JSON"),
    ])
    def test_refusals(self, staged, choices, logged_in, message):
        session_id, share = staged
        body = self.run(session_id, choices, logged_in=logged_in).json()
        assert body["success"] == 0 and body["message"].startswith(message)
        assert json.loads((share / "metadata.json").read_text())["questionable_obs_columns"] == QUESTIONABLE

    def test_no_metadata(self, staged):
        session_id, _ = staged
        body = self.run(session_id, {"replicate": "categorical"}, share_uid="other").json()
        assert body == {"success": 0, "message": "No metadata JSON file found for this dataset."}

    @pytest.mark.parametrize("dataset_format, choices, fragment", [
        ("single-cell-rnaseq", {"replicate": "categorical"}, "Dataset file not found"),
        ("spatial", {"replicate": "categorical"}, "Spatial dataset file not found"),
        ("single-cell-rnaseq", {"replicate": "ordinal"}, "Unknown dtype choice"),
    ])
    def test_failures_leave_dataset_unreviewed(self, staged, dataset_format, choices, fragment):
        session_id, share = staged
        if choices["replicate"] == "ordinal":
            stage_h5ad(share)
        metadata = {"dataset_format": dataset_format, "questionable_obs_columns": QUESTIONABLE}
        (share / "metadata.json").write_text(json.dumps(metadata))
        body = self.run(session_id, choices).json()
        assert body["success"] == 0 and body["message"].startswith("Error applying dtype choices:")
        assert fragment in body["message"]
        assert json.loads((share / "metadata.json").read_text()) == metadata


# ---------------------------------------------------------------- dataset obs listings

@pytest.fixture
def dataset_h5ad(tmp_path):
    """Write an h5ad with the given obs columns (4 cells) and return (db spec, path)."""
    anndata = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    import numpy as np
    import pandas as pd

    def make(columns, n_cells=4):
        path = tmp_path / "DS1.h5ad"
        obs = pd.DataFrame(columns, index=[f"c{i}" for i in range(n_cells)])
        anndata.AnnData(X=np.ones((n_cells, 2), dtype=np.float32), obs=obs).write_h5ad(path)
        return {"datasets": [{"id": "DS1", "h5ad_path": str(path)}]}
    return make


class TestGetH5adObsColumns:
    def test_lists_columns_without_replicate(self, dataset_h5ad):
        db = dataset_h5ad({"cell_type": ["a", "b", "a", "b"], "replicate": ["1", "1", "2", "2"], "age": [1, 2, 3, 4]})
        assert run_cgi("get_h5ad_obs_columns.cgi", db=db, query={"dataset_id": "DS1"}).json() == \
            {"success": 1, "obs_columns": ["cell_type", "age"]}

    def test_missing_file(self):
        result = run_cgi("get_h5ad_obs_columns.cgi", db={"datasets": [{"id": "DS1"}]}, query={"dataset_id": "DS1"})
        assert result.json() == {"success": 0, "error": "No file found for this dataset"}


class TestGetConditionList:
    def run(self, db):
        return run_cgi("get_condition_list.cgi", db=db, query={"dataset_id": "DS1"}).json()

    def test_conditions_from_obs_columns(self, dataset_h5ad):
        # replicate is reported, not combined; louvain and barcode are skipped
        db = dataset_h5ad({"cell_type": ["a", "b", "a", "b"], "age": ["P1", "P1", "P2", "P1"],
                           "replicate": ["1", "1", "2", "2"], "louvain": ["0", "1", "0", "1"],
                           "barcode": ["w", "x", "y", "z"]})
        assert self.run(db) == {"has_replicates": 1, "conditions": [
            {"class_label": "a;P1", "formatted_class_label": "a, P1"},
            {"class_label": "b;P1", "formatted_class_label": "b, P1"},
            {"class_label": "a;P2", "formatted_class_label": "a, P2"}]}

    def test_too_many_columns(self, dataset_h5ad):
        db = dataset_h5ad({f"col{i}": ["a", "b", "a", "b"] for i in range(6)})
        body = self.run(db)
        assert body["success"] == 0 and "too many observation columns" in body["error"]

    def test_no_columns_on_large_dataset(self, dataset_h5ad):
        assert self.run(dataset_h5ad({}, n_cells=51)) == {"has_replicates": 0, "success": -1}

    def test_missing_file(self):
        assert self.run({"datasets": [{"id": "DS1"}]}) == {"has_replicates": 0, "success": 0,
                                                            "error": "No file found for DS1."}
