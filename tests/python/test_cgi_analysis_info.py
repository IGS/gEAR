"""Workbench analysis lookups: get_embedded_tsne_display.cgi responses and Analysis.discover_vetting()."""

from pathlib import Path

import pytest

import gear.analysis as analysis
from helpers.cgi_harness import run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
WORKBENCH_ANALYSIS = {"gear.analysis": FAKES / "workbench_analysis.py"}
DB = {"datasets": [{"id": "DS1", "owner_id": 1}]}


class TestEmbeddedTsneDisplay:
    def run(self, tmp_path, **query):
        return run_cgi("get_embedded_tsne_display.cgi", db=DB, query=query, fake_modules=WORKBENCH_ANALYSIS,
                       env={"GEAR_TEST_ANALYSIS_DIR": str(tmp_path)})

    def test_returns_config_from_primary_analysis(self, tmp_path):
        pytest.importorskip("anndata")
        body = self.run(tmp_path, dataset_id="DS1").json()
        # The fake dataset has a subclass_label column but no tSNE coordinates
        assert body == {"plotly_config": {"colorize_legend_by": "subclass_label", "x_axis": None, "y_axis": None}}

    @pytest.mark.parametrize("query, error", [({}, "No dataset_id provided."),
                                              ({"dataset_id": "NOPE"}, "No dataset found with that ID.")])
    def test_errors_are_returned_as_json(self, tmp_path, query, error):
        # These used to be printed to /dev/null, giving "End of script output before headers"
        assert self.run(tmp_path, **query).json() == {"error": error}


class TestDiscoverVetting:
    def make(self, ana_type, user_id):
        return analysis.Analysis(id="A1", dataset_id="DS1", type=ana_type, user_id=user_id)

    @pytest.fixture(autouse=True)
    def users(self, monkeypatch):
        users = {1: type("User", (), {"is_curator": 0})(), 2: type("User", (), {"is_curator": 1})()}
        monkeypatch.setattr(analysis, "get_user_by_id", lambda user_id: users.get(user_id))

    def test_primary_is_gear_vetted(self):
        assert self.make("primary", None).discover_vetting(current_user_id=1) == "gear"

    def test_unknown_owner_leaves_vetting_unset(self):
        assert self.make("public", None).discover_vetting(current_user_id=1) is None

    @pytest.mark.parametrize("current_user_id, expected", [(1, "owner"), (2, "gear"), (3, None)])
    def test_owned_analysis(self, current_user_id, expected):
        ana = self.make("user_saved", 1)
        if expected is None:
            with pytest.raises(Exception, match="without a current user"):
                ana.discover_vetting(current_user_id=current_user_id)
        else:
            assert ana.discover_vetting(current_user_id=current_user_id) == expected
