"""gear.dominosignal: finding DominoSignal outputs, and reading the dataset details from each output's matching .json file."""

import json

import pytest
from gear import dominosignal as ds

FULL_INFO = {
    "Title": " Intermittent clearance of p21-highly-expressing cells extends lifespan ",
    "Author": "B Wang",
    "Year": "2024",
    "DOI": "10.1016/j.cmet.2024.07.006",
    "PMID": "39111286",
    "gEAR Dataset URL": "umgear.org",
    "Description": "P21 Dataset",
}


@pytest.fixture
def domino_dir(tmp_path, monkeypatch):
    """An empty DominoSignal output directory, in place of www/cache/domino_signal."""
    monkeypatch.setattr(ds, "DOMINO_DIR", str(tmp_path))
    return tmp_path


def add_output(domino_dir, name, info=None, raw_json=None):
    """Adds <name>.parquet (empty; only its existence matters here), with a <name>.json of info or raw_json if given."""
    (domino_dir / f"{name}.parquet").touch()
    if info is not None:
        raw_json = json.dumps(info)
    if raw_json is not None:
        (domino_dir / f"{name}.json").write_text(raw_json)
    return f"{name}.parquet"


def test_lists_only_parquet_files(domino_dir):
    add_output(domino_dir, "b_output", FULL_INFO)
    add_output(domino_dir, "a_output")
    (domino_dir / "notes.txt").touch()
    assert ds.list_domino_files() == ["a_output.parquet", "b_output.parquet"]


def test_reads_dataset_info_from_matching_json(domino_dir):
    filename = add_output(domino_dir, "domino_network_p21.reformatted", FULL_INFO)
    assert ds.load_dataset_info(filename) == {
        "title": "Intermittent clearance of p21-highly-expressing cells extends lifespan",  # whitespace trimmed
        "author": "B Wang",
        "year": "2024",
        "doi": "10.1016/j.cmet.2024.07.006",
        "pmid": "39111286",
        "gear_url": "umgear.org",
        "description": "P21 Dataset",
    }


def test_blank_and_missing_fields_are_none(domino_dir):
    filename = add_output(domino_dir, "partial", {"Title": "A paper", "Author": "   ", "Year": "", "DOI": None})
    info = ds.load_dataset_info(filename)
    assert info["title"] == "A paper"
    assert info["author"] is None and info["year"] is None and info["doi"] is None
    assert info["pmid"] is None and info["gear_url"] is None and info["description"] is None


def test_keys_are_case_insensitive_and_numbers_become_strings(domino_dir):
    filename = add_output(domino_dir, "lowercase", {"title": "A paper", "AUTHOR": "X Y", "year": 2023, "pmid": 123,
                                                    "gear dataset url": "https://umgear.org/p?share_id=1"})
    info = ds.load_dataset_info(filename)
    assert (info["title"], info["author"], info["year"], info["pmid"]) == ("A paper", "X Y", "2023", "123")
    assert info["gear_url"] == "https://umgear.org/p?share_id=1"


@pytest.mark.parametrize("raw_json", [None, "{ not json", "[1, 2]", '"just a string"'])
def test_missing_or_unusable_json_gives_none(domino_dir, raw_json):
    filename = add_output(domino_dir, "no_info", raw_json=raw_json)
    assert ds.load_dataset_info(filename) is None


def test_dataset_info_only_reads_json_beside_a_known_output(domino_dir):
    # The filename is checked exactly as for the parquet itself, so it can't reach outside the directory
    with pytest.raises(ds.DominoFileError):
        ds.load_dataset_info("../outside.parquet")
    with pytest.raises(ds.DominoFileError):
        ds.load_dataset_info("missing.parquet")
    with pytest.raises(ds.DominoFileError):
        ds.load_dataset_info("domino.json")
