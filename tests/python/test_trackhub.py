"""lib/gear/trackhub.py: hub/genomes/trackDb parsing, hub validation and the public-URL guard."""

import importlib
import socket
import sys
import types
from pathlib import Path

import pytest

REAL_TRACKHUB = Path(__file__).resolve().parents[2] / "lib" / "gear" / "trackhub.py"

HUB_JSON = {"hub": "myhub", "shortLabel": "My hub", "longLabel": "My hub", "email": "a@b.org",
            "useOneFile": "on", "genome": "mm10"}


def _track(**fields):
    return {"track": "t1", "type": "bigWig", "shortLabel": "T1", "longLabel": "Track 1",
            "visibility": "full", "bigDataUrl": "https://example.org/t1.bw", **fields}


@pytest.fixture(scope="module")
def trackhub():
    """The real gear.trackhub, even if another test left the stub from fakes/api_trackhub.py in sys.modules."""
    pytest.importorskip("pyBigWig")
    pytest.importorskip("Bio")
    with pytest.MonkeyPatch.context() as mp:
        module = sys.modules.get("gear.trackhub")
        if module is None or Path(getattr(module, "__file__", "")).resolve() != REAL_TRACKHUB:
            mp.delitem(sys.modules, "gear.trackhub", raising=False)
            module = importlib.import_module("gear.trackhub")
        assert Path(module.__file__).resolve() == REAL_TRACKHUB
        yield module


@pytest.fixture
def addrinfo(trackhub, monkeypatch):
    """Resolve every host to the given IPs, recording the lookups."""
    lookups = []

    def install(*ips):
        def getaddrinfo(host, port, *args, **kwargs):
            lookups.append(host)
            return [(socket.AF_INET6 if ":" in ip else socket.AF_INET, socket.SOCK_STREAM, 6, "", (ip, 0)) for ip in ips]
        monkeypatch.setattr(trackhub.socket, "getaddrinfo", getaddrinfo)
        return lookups

    monkeypatch.delenv("ENVIRONMENT", raising=False)
    return install


@pytest.fixture
def no_http(trackhub, monkeypatch):
    """Fail on any HEAD request; tests that expect one replace it."""
    def head(*args, **kwargs):
        raise AssertionError("unexpected requests.head")
    monkeypatch.setattr(trackhub.requests, "head", head)


# ---------------------------------------------------------------- fetch_trackdb_path

GENOMES_TXT = """\
genome mm10
trackDb mm10/trackDb.txt

genome hg38
trackDb hg38/trackDb.txt
"""


@pytest.mark.parametrize("assembly, expected", [("mm10", "mm10/trackDb.txt"), ("hg38", "hg38/trackDb.txt")])
def test_fetch_trackdb_path(trackhub, assembly, expected):
    assert trackhub.fetch_trackdb_path(GENOMES_TXT, assembly) == expected


@pytest.mark.parametrize("assembly", ["mm1", "mm", "hg3", "mm39"])
def test_fetch_trackdb_path_matches_whole_assembly(trackhub, assembly):
    with pytest.raises(ValueError, match="not found"):
        trackhub.fetch_trackdb_path(GENOMES_TXT, assembly)


def test_fetch_trackdb_path_prefix_assembly_finds_its_own_stanza(trackhub):
    genomes = "genome mm10\ntrackDb mm10/trackDb.txt\n\ngenome mm1\ntrackDb mm1/trackDb.txt\n"
    assert trackhub.fetch_trackdb_path(genomes, "mm1") == "mm1/trackDb.txt"


def test_fetch_trackdb_path_trackdb_anywhere_in_stanza(trackhub):
    genomes = """\
genome hg38
twoBitPath hg38/hg38.2bit
description Human
trackDb hg38/trackDb.txt

genome mm10
trackDb mm10/trackDb.txt
"""
    assert trackhub.fetch_trackdb_path(genomes, "hg38") == "hg38/trackDb.txt"


def test_fetch_trackdb_path_does_not_cross_stanzas(trackhub):
    genomes = "genome hg38\ndescription Human\n\ntrackDb other/trackDb.txt\n"
    with pytest.raises(ValueError):
        trackhub.fetch_trackdb_path(genomes, "hg38")


def test_fetch_trackdb_path_missing_trackdb(trackhub):
    with pytest.raises(ValueError, match="mm10"):
        trackhub.fetch_trackdb_path("genome mm10\n", "mm10")


# ---------------------------------------------------------------- _is_safe_public_http_url

def test_is_safe_public_http_url_public_ip(trackhub, addrinfo):
    lookups = addrinfo("93.184.216.34", "2606:2800:220:1:248:1893:25c8:1946")
    assert trackhub._is_safe_public_http_url("https://example.org/a.bw") is True
    assert lookups == ["example.org"]


@pytest.mark.parametrize("ip", ["127.0.0.1", "10.1.2.3", "172.16.0.1", "192.168.0.10", "169.254.169.254",
                                "0.0.0.0", "224.0.0.1", "240.0.0.1", "::1", "fe80::1", "fd00::1"])
def test_is_safe_public_http_url_rejects_non_public_ip(trackhub, addrinfo, ip):
    addrinfo(ip)
    assert trackhub._is_safe_public_http_url("https://example.org/a.bw") is False


def test_is_safe_public_http_url_rejects_if_any_address_is_private(trackhub, addrinfo):
    addrinfo("93.184.216.34", "10.0.0.1")
    assert trackhub._is_safe_public_http_url("https://example.org/a.bw") is False


@pytest.mark.parametrize("url", ["ftp://example.org/a.bw", "file:///etc/passwd", "https:///a.bw", "not a url"])
def test_is_safe_public_http_url_rejects_bad_urls_without_lookup(trackhub, addrinfo, url):
    lookups = addrinfo("93.184.216.34")
    assert trackhub._is_safe_public_http_url(url) is False
    assert lookups == []


def test_is_safe_public_http_url_lookup_failure(trackhub, monkeypatch):
    def getaddrinfo(*args, **kwargs):
        raise socket.gaierror("no such host")
    monkeypatch.setattr(trackhub.socket, "getaddrinfo", getaddrinfo)
    monkeypatch.delenv("ENVIRONMENT", raising=False)
    assert trackhub._is_safe_public_http_url("https://nowhere.invalid/a.bw") is False


def test_is_safe_public_http_url_development_skips_lookup(trackhub, addrinfo, monkeypatch):
    lookups = addrinfo("127.0.0.1")
    monkeypatch.setenv("ENVIRONMENT", "Development")
    assert trackhub._is_safe_public_http_url("http://localhost:8080/a.bw") is True
    assert trackhub._is_safe_public_http_url("ftp://localhost/a.bw") is False
    assert lookups == []


# ---------------------------------------------------------------- _normalize_track_dict

def test_normalize_track_dict_keeps_known_and_gos_fields(trackhub, no_http):
    raw = {"track": "t1", "type": "bigWig", "color": "31,119,180", "gos_mark": "bar", "html": "x.html",
           "bigDataUrl": "t1.bw"}
    assert trackhub._normalize_track_dict(raw) == {
        "track": "t1", "type": "bigWig", "color": "rgb(31,119,180)", "gos_mark": "bar", "bigDataUrl": "t1.bw",
    }


def test_normalize_track_dict_keeps_rgb_color(trackhub):
    assert trackhub._normalize_track_dict({"track": "t1", "color": "rgb(1,2,3)"})["color"] == "rgb(1,2,3)"


def test_normalize_track_dict_requires_track(trackhub):
    with pytest.raises(trackhub.TrackHubProcessingError, match="'track' field"):
        trackhub._normalize_track_dict({"type": "bigWig"})


def test_normalize_track_dict_container_skips_url_resolution(trackhub, no_http):
    raw = {"track": "c1", "container": "multiWig", "bigDataUrl": "c1.bw"}
    assert trackhub._normalize_track_dict(raw, resolve_urls=True, trackdb_url="https://example.org/") == raw


def test_normalize_track_dict_resolves_relative_url_and_follows_redirect(trackhub, addrinfo, monkeypatch):
    addrinfo("93.184.216.34")
    heads = []

    def head(url, **kwargs):
        heads.append((url, kwargs))
        return types.SimpleNamespace(status_code=200, url="https://cdn.example.org/t1.bw")
    monkeypatch.setattr(trackhub.requests, "head", head)

    track = trackhub._normalize_track_dict({"track": "t1", "bigDataUrl": "data/t1.bw"}, resolve_urls=True,
                                           trackdb_url="https://example.org/hub/mm10/trackDb.txt")
    assert heads == [("https://example.org/hub/mm10/data/t1.bw", {"allow_redirects": True, "timeout": 10})]
    assert track["bigDataUrl"] == "https://cdn.example.org/t1.bw"


def test_normalize_track_dict_ignores_redirect_to_private_host(trackhub, monkeypatch):
    monkeypatch.delenv("ENVIRONMENT", raising=False)
    monkeypatch.setattr(trackhub.socket, "getaddrinfo",
                        lambda host, port: [(socket.AF_INET, socket.SOCK_STREAM, 6, "",
                                             ("10.0.0.1" if host == "internal" else "93.184.216.34", 0))])
    monkeypatch.setattr(trackhub.requests, "head",
                        lambda url, **kwargs: types.SimpleNamespace(status_code=200, url="http://internal/t1.bw"))
    track = trackhub._normalize_track_dict({"track": "t1", "bigDataUrl": "https://example.org/t1.bw"},
                                           resolve_urls=True, trackdb_url="https://example.org/trackDb.txt")
    assert track["bigDataUrl"] == "https://example.org/t1.bw"


def test_normalize_track_dict_does_not_request_unsafe_url(trackhub, addrinfo, no_http):
    addrinfo("169.254.169.254")
    track = trackhub._normalize_track_dict({"track": "t1", "bigDataUrl": "http://metadata/t1.bw"},
                                           resolve_urls=True, trackdb_url="https://example.org/trackDb.txt")
    assert track["bigDataUrl"] == "http://metadata/t1.bw"


@pytest.mark.parametrize("head_result", [
    types.SimpleNamespace(status_code=404, url="https://example.org/missing.bw"),
    OSError("connection refused"),
])
def test_normalize_track_dict_keeps_url_when_head_fails(trackhub, addrinfo, monkeypatch, head_result):
    addrinfo("93.184.216.34")

    def head(url, **kwargs):
        if isinstance(head_result, Exception):
            raise head_result
        return head_result
    monkeypatch.setattr(trackhub.requests, "head", head)
    track = trackhub._normalize_track_dict({"track": "t1", "bigDataUrl": "t1.bw"}, resolve_urls=True,
                                           trackdb_url="https://example.org/hub/trackDb.txt")
    assert track["bigDataUrl"] == "https://example.org/hub/t1.bw"


# ---------------------------------------------------------------- stanza parsing

TRACKDB_TXT = """\
# comment
track multi
container multiWig
type bigWig
shortLabel Multi

track t1
parent multi
type bigWig
bigDataUrl https://example.org/t1.bw
shortLabel Track 1
longLabel Track one, replicate 1
color 255,0,0
visibility dense
gos_height 80
priority 3
"""


def test_parse_track_stanzas(trackhub, no_http):
    tracks = trackhub._parse_track_stanzas(TRACKDB_TXT.splitlines())
    assert tracks == [
        {"track": "multi", "container": "multiWig", "type": "bigWig", "shortLabel": "Multi"},
        {"track": "t1", "parent": "multi", "type": "bigWig", "bigDataUrl": "https://example.org/t1.bw",
         "shortLabel": "Track 1", "longLabel": "Track one, replicate 1", "color": "rgb(255,0,0)",
         "visibility": "dense", "gos_height": "80"},
    ]


def test_parse_track_stanzas_ignores_fields_before_first_track(trackhub):
    assert trackhub._parse_track_stanzas(["type bigWig", "gos_height 10", "", "track t1"]) == [{"track": "t1"}]


def test_parse_track_stanzas_empty(trackhub):
    assert trackhub._parse_track_stanzas([]) == []


def test_parse_tracks_from_trackdb_resolves_relative_urls(trackhub, monkeypatch):
    monkeypatch.setenv("ENVIRONMENT", "development")
    monkeypatch.setattr(trackhub.requests, "head",
                        lambda url, **kwargs: types.SimpleNamespace(status_code=200, url=url))
    tracks = trackhub.parse_tracks_from_trackdb("track t1\ntype bigWig\nbigDataUrl ../data/t1.bw\n",
                                                "http://localhost:8080/hubs/h1/mm10/trackDb.txt")
    assert tracks == [{"track": "t1", "type": "bigWig", "bigDataUrl": "http://localhost:8080/hubs/h1/data/t1.bw"}]


HUB_TXT = """\
hub myhub
shortLabel My hub
longLabel My hub, longer
useOneFile on
email someone@example.org
descriptionUrl about.html

genome mm10

track t1
type bigWig
bigDataUrl t1.bw
shortLabel T1

track t2
type bigBed
bigDataUrl https://example.org/t2.bb
"""


def test_parse_hub_from_file(trackhub, no_http, tmp_path):
    hub_file = tmp_path / "hub.txt"
    hub_file.write_text(HUB_TXT)
    hub_json, tracks = trackhub.parse_hub_from_file(hub_file.read_text())
    assert hub_json == {"hub": "myhub", "shortLabel": "My hub", "longLabel": "My hub, longer", "useOneFile": "on",
                        "email": "someone@example.org", "genome": "mm10"}
    assert tracks == [
        {"track": "t1", "type": "bigWig", "bigDataUrl": "t1.bw", "shortLabel": "T1"},
        {"track": "t2", "type": "bigBed", "bigDataUrl": "https://example.org/t2.bb"},
    ]


def test_parse_hub_from_file_track_labels_do_not_leak_into_hub(trackhub):
    hub_json, tracks = trackhub.parse_hub_from_file("hub h\nshortLabel Hub\n\ntrack t1\nshortLabel Track\n")
    assert hub_json == {"hub": "h", "shortLabel": "Hub"}
    assert tracks == [{"track": "t1", "shortLabel": "Track"}]


def test_parse_hub_from_file_without_tracks(trackhub):
    assert trackhub.parse_hub_from_file("hub h\ngenomesFile genomes.txt\n") == ({"hub": "h", "genomesFile": "genomes.txt"}, [])


# ---------------------------------------------------------------- validate_hub_contents

def test_validate_hub_contents_accepts_valid_hub(trackhub):
    container = _track(track="c1", container="multiWig")
    del container["bigDataUrl"]
    trackhub.validate_hub_contents(HUB_JSON, [container, _track(parent="c1"), _track(track="t2", type="hic")])


@pytest.mark.parametrize("field", ["hub", "shortLabel", "longLabel", "email", "useOneFile", "genome"])
def test_validate_hub_contents_missing_hub_field(trackhub, field):
    hub_json = {k: v for k, v in HUB_JSON.items() if k != field}
    with pytest.raises(trackhub.TrackHubProcessingError, match=f"'{field}' field"):
        trackhub.validate_hub_contents(hub_json, [])


@pytest.mark.parametrize("field", ["track", "type", "shortLabel", "longLabel", "visibility"])
def test_validate_hub_contents_missing_track_field(trackhub, field):
    track = _track()
    del track[field]
    with pytest.raises(trackhub.TrackHubProcessingError, match=f"'{field}' field"):
        trackhub.validate_hub_contents(HUB_JSON, [track])


@pytest.mark.parametrize("tracks, message", [
    ([{k: v for k, v in _track().items() if k != "bigDataUrl"}], "either a 'bigDataUrl'"),
    ([_track(type="bam")], "unsupported type 'bam'"),
    ([_track(container="composite")], "unsupported container type 'composite'"),
    ([_track(parent="nope")], "parent track 'nope' that doesn't exist"),
    ([_track(track="p"), _track(track="t2", parent="p")], "isn't a valid container"),
])
def test_validate_hub_contents_rejects(trackhub, tracks, message):
    with pytest.raises(trackhub.TrackHubProcessingError, match=message):
        trackhub.validate_hub_contents(HUB_JSON, tracks)


def test_validate_hub_contents_drops_null_container_and_parent(trackhub):
    track = _track(container=None, parent=None)
    trackhub.validate_hub_contents(HUB_JSON, [track])
    assert "container" not in track and "parent" not in track


@pytest.mark.parametrize("genome", ["hg19", "HG38"])
@pytest.mark.parametrize("track_type", ["vcfTabix", "hic"])
def test_validate_hub_contents_human_privacy_types(trackhub, genome, track_type):
    with pytest.raises(trackhub.TrackHubProcessingError, match="human genome"):
        trackhub.validate_hub_contents({**HUB_JSON, "genome": genome}, [_track(type=track_type)])


def test_validate_hub_contents_human_allows_bigwig(trackhub):
    trackhub.validate_hub_contents({**HUB_JSON, "genome": "hg38"}, [_track()])
