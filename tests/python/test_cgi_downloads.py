"""Download CGIs: "Is downloadable" enforcement, HTTP status codes and attachment file names."""

from pathlib import Path

import pytest

from helpers.cgi_harness import run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
ANALYSIS = {"gear.analysis": FAKES / "downloads_analysis.py"}

# Session "owner" is user 1, who owns both datasets; "other" is user 2.
SESSIONS = {"owner": 1, "other": 2}
PRIVATE = {"id": "D1", "share_id": "priv", "owner_id": 1, "is_downloadable": 0}
OPEN = {"id": "D2", "share_id": "open", "owner_id": 1, "is_downloadable": 1}


def db_with(*datasets, **extra):
    return {"sessions": SESSIONS, "datasets": list(datasets), **extra}


class TestDownloadSourceFile:
    def run(self, db=None, cookies=None, **query):
        return run_cgi("download_source_file.cgi", query=query, cookies=cookies, fake_modules=ANALYSIS,
                       db=db if db is not None else db_with(PRIVATE, OPEN))

    def test_missing_type(self):
        result = self.run(share_id="open")
        assert result.status == 400
        assert result.headers["content-type"] == "text/plain"
        assert "Type must be provided" in result.text

    def test_missing_ids(self):
        result = self.run(type="tarball")
        assert result.status == 400

    def test_unknown_share(self):
        result = self.run(share_id="nope", type="tarball")
        assert result.status == 404
        assert "nope" in result.text

    @pytest.mark.parametrize("cookies", [None, {"gear_session_id": "other"}])
    def test_non_downloadable_refused_for_non_owner(self, cookies):
        result = self.run(share_id="priv", type="tarball", cookies=cookies)
        assert result.status == 403
        assert "not available for download" in result.text

    def test_non_downloadable_refused_for_other_user_session_param(self):
        result = self.run(share_id="priv", type="h5ad", session_id="other")
        assert result.status == 403

    @pytest.mark.parametrize("how", ["cookie", "param"])
    def test_owner_passes_the_check(self, how):
        if how == "cookie":
            result = self.run(share_id="priv", type="tarball", cookies={"gear_session_id": "owner"})
        else:
            result = self.run(share_id="priv", type="tarball", session_id="owner")
        assert result.status == 404
        assert "File not found" in result.text

    def test_downloadable_passes_for_anonymous(self):
        result = self.run(share_id="open", type="h5ad")
        assert result.status == 404
        assert "File not found" in result.text

    @pytest.mark.parametrize("dataset, cookies", [
        (PRIVATE, {"gear_session_id": "owner"}),
        (OPEN, None),
    ])
    def test_tarball_download(self, tmp_path, dataset, cookies):
        archive = tmp_path / "x.tar.gz"
        archive.write_bytes(b"\x1f\x8b fake tarball bytes")
        result = self.run(db=db_with({**dataset, "archive_path": str(archive)}),
                          share_id=dataset["share_id"], type="tarball", cookies=cookies)
        assert result.status == 200
        assert result.headers["content-disposition"] == f"attachment; filename={dataset['share_id']}.tar.gz"
        assert result.body == archive.read_bytes()

    def test_h5ad_download_for_owner(self, tmp_path):
        h5ad = tmp_path / "D1.h5ad"
        h5ad.write_bytes(b"HDF5 fake")
        result = self.run(db=db_with({**PRIVATE, "h5ad_path": str(h5ad)}),
                          share_id="priv", type="h5ad", session_id="owner")
        assert result.status == 200
        assert result.headers["content-disposition"] == "attachment; filename=priv.h5ad"
        assert result.body == b"HDF5 fake"

    def test_h5ad_refused_for_non_owner_even_if_file_exists(self, tmp_path):
        h5ad = tmp_path / "D1.h5ad"
        h5ad.write_bytes(b"HDF5 fake")
        result = self.run(db=db_with({**PRIVATE, "h5ad_path": str(h5ad)}), share_id="priv", type="h5ad")
        assert result.status == 403
        assert b"HDF5" not in result.body

    def test_metadata_not_blocked_for_non_downloadable(self):
        result = self.run(db=db_with(PRIVATE, metadata={"priv": "a,b\n1,2\n"}), share_id="priv", type="metadata")
        assert result.status == 200
        assert result.headers["content-disposition"] == "attachment; filename=priv.metadata.csv"
        assert result.text == "a,b\n1,2\n"

    def test_missing_metadata(self):
        result = self.run(share_id="priv", type="metadata")
        assert result.status == 404
        assert "Metadata not found" in result.text


class TestDownloadProjection:
    def run(self, **query):
        return run_cgi("download_projection.cgi", query=query, db=db_with(OPEN))

    def test_no_ids(self):
        result = self.run(projection_id="P1")
        assert result.status == 400
        assert result.headers["content-type"] == "text/plain"

    def test_unknown_share(self):
        result = self.run(share_id="nope", projection_id="P1")
        assert result.status == 404

    def test_no_output_file(self):
        result = self.run(share_id="open", projection_id="does-not-exist")
        assert result.status == 404
        assert "Projection output not found" in result.text


class TestDownloadWeightedGeneCart:
    def test_missing_share_id(self):
        result = run_cgi("download_weighted_gene_cart.cgi")
        assert result.status == 400

    def test_unknown_cart(self):
        result = run_cgi("download_weighted_gene_cart.cgi", query={"share_id": "pytest-does-not-exist"})
        assert result.status == 404

    def test_existing_cart(self, tmp_cart):
        content = "id\tgene\tP1\nENSG1\tSox2\t0.5\n"
        share_id = tmp_cart(content)
        result = run_cgi("download_weighted_gene_cart.cgi", query={"share_id": share_id})
        assert result.status == 200
        assert result.headers["content-disposition"] == f"attachment; filename=cart.{share_id}.tab"
        assert result.body == content.encode()


class TestGetPatternWeightedGenes:
    CART = "id\tgene_symbol\tP1\tP2\nE1\tSox2\t0.1\t9\nE2\tPax6\t0.9\t8\nE3\tGata3\t0.5\t7\n"

    def test_path_traversal_rejected(self):
        result = run_cgi("get_pattern_weighted_genes.cgi", query={"source_id": "../../x", "pattern_id": "P1"})
        assert result.status == 400
        body = result.json()
        assert body["success"] == 0 and body["error"]

    def test_unknown_cart(self):
        result = run_cgi("get_pattern_weighted_genes.cgi",
                         query={"source_id": "pytest-does-not-exist", "pattern_id": "P1"})
        assert result.status == 404
        assert result.json()["success"] == 0

    def test_unknown_pattern(self, tmp_cart):
        share_id = tmp_cart(self.CART)
        result = run_cgi("get_pattern_weighted_genes.cgi", query={"source_id": share_id, "pattern_id": "P9"})
        assert result.status == 404
        assert "P9" in result.json()["error"]

    def test_found_sorted_by_weight(self, tmp_cart):
        share_id = tmp_cart(self.CART)
        result = run_cgi("get_pattern_weighted_genes.cgi", query={"source_id": share_id, "pattern_id": "P1"})
        assert result.status == 200
        assert result.json() == [
            {"gene": "Pax6", "weight": 0.9},
            {"gene": "Gata3", "weight": 0.5},
            {"gene": "Sox2", "weight": 0.1},
        ]
