"""resources.gosling_spec helpers: hub URL validation (SSRF guard), position strings, assembly arrays and URL checks."""

import types

import pytest

pytest.importorskip("gosling")
from resources import gosling_spec  # noqa: E402

DOMAIN_URL = "https://gear.example.org"


# ---------------------------------------------------------------- _validate_hub_url

@pytest.fixture
def domain(monkeypatch):
    monkeypatch.setattr(gosling_spec.geardb, "_read_domain_url", lambda: DOMAIN_URL)
    monkeypatch.delenv("ENVIRONMENT", raising=False)


@pytest.fixture
def development(domain, monkeypatch):
    monkeypatch.setenv("ENVIRONMENT", "development")


@pytest.mark.parametrize("url", [
    "https://gear.example.org/hub.txt",
    "http://GEAR.example.org/hub.txt",
    "https://tracks.gear.example.org/hubs/hub.txt",
])
def test_validate_hub_url_accepts(domain, url):
    gosling_spec._validate_hub_url(url)


@pytest.mark.parametrize("url, message", [
    ("ftp://gear.example.org/hub.txt", "scheme"),
    ("file:///etc/passwd", "scheme"),
    ("https:///hub.txt", "missing hostname"),
    ("http://127.0.0.1/hub.txt", "private"),
    ("http://10.0.0.5/hub.txt", "private"),
    ("http://192.168.1.1/hub.txt", "private"),
    ("http://169.254.169.254/latest/meta-data", "private"),
    ("http://[::1]/hub.txt", "private"),
    ("https://example.org/hub.txt", "not in the allowed list"),
    ("https://evilgear.example.org/hub.txt", "not in the allowed list"),
    ("https://gear.example.org.evil.com/hub.txt", "not in the allowed list"),
    ("https://umgear.org/hub.txt", "not in the allowed list"),
    ("http://localhost/hub.txt", "not in the allowed list"),
])
def test_validate_hub_url_rejects(domain, url, message):
    with pytest.raises(ValueError, match=message):
        gosling_spec._validate_hub_url(url)


@pytest.mark.parametrize("url", ["https://8.8.8.8/hub.txt", "http://[2001:4860:4860::8888]/hub.txt"])
def test_validate_hub_url_rejects_public_ip_outside_allowlist(domain, url):
    with pytest.raises(ValueError, match="not in the allowed list"):
        gosling_spec._validate_hub_url(url)


@pytest.mark.parametrize("domain_url", ["https://203.0.113.7", "http://10.1.2.3:8080"])
def test_validate_hub_url_accepts_configured_ip(monkeypatch, domain_url):
    """A server with no DNS name yet is configured by IP address, and its own hubs must still load."""
    monkeypatch.setattr(gosling_spec.geardb, "_read_domain_url", lambda: domain_url)
    monkeypatch.delenv("ENVIRONMENT", raising=False)
    host = domain_url.split("//")[1]
    gosling_spec._validate_hub_url(f"http://{host}/hub.txt")
    with pytest.raises(ValueError):
        gosling_spec._validate_hub_url("https://203.0.113.8/hub.txt")


def test_validate_hub_url_rejects_domain_url_without_scheme(monkeypatch):
    monkeypatch.setattr(gosling_spec.geardb, "_read_domain_url", lambda: "gear.example.org")
    with pytest.raises(ValueError, match="not a valid URL"):
        gosling_spec._validate_hub_url("https://gear.example.org/hub.txt")


def test_validate_hub_url_requires_configured_domain(monkeypatch):
    monkeypatch.setattr(gosling_spec.geardb, "_read_domain_url", lambda: None)
    with pytest.raises(ValueError, match="not configured"):
        gosling_spec._validate_hub_url("https://gear.example.org/hub.txt")


@pytest.mark.parametrize("url", [
    "https://umgear.org/hub.txt",
    "https://devel.umgear.org/hub.txt",
    "http://web/hub.txt",
    "http://localhost:8080/hub.txt",
    "http://127.0.0.1:8080/hub.txt",
    "https://gear.example.org/hub.txt",
])
def test_validate_hub_url_development_extras(development, url):
    gosling_spec._validate_hub_url(url)


@pytest.mark.parametrize("url", [
    # Each ends with a single character of "umgear.org"
    "https://example.org/hub.txt",
    "https://example.com/hub.txt",
    "https://evilumgear.org/hub.txt",
    "https://notweb/hub.txt",
    "http://127.0.0.2/hub.txt",
    "http://0.0.0.0/hub.txt",
    "http://169.254.169.254/latest/meta-data",
    "https://8.8.8.8/hub.txt",
])
def test_validate_hub_url_development_still_restricted(development, url):
    with pytest.raises(ValueError):
        gosling_spec._validate_hub_url(url)


# ---------------------------------------------------------------- position strings

def test_parse_position_str():
    assert gosling_spec.parse_position_str("mm10.chr7:142575529-142578143") == ("mm10", "chr7", 142575529, 142578143)


# Dotted (non-UCSC) assembly names and comma-grouped coordinates are never produced internally, so are unsupported
@pytest.mark.parametrize("position", [
    "NA", "", "mm10.chr7", "mm10.chr7:100", "mm10.chr7:a-b", "GRCh38.p14.chr1:100-200", "hg38.chr1:1,000-2,000",
])
def test_parse_position_str_invalid_returns_nones(position):
    assert gosling_spec.parse_position_str(position) == (None, None, None, None)


class _View:
    def properties(self, **kwargs):
        return types.SimpleNamespace(**kwargs)


@pytest.mark.parametrize("hic_found, padding", [(False, 1500), (True, 1500 * 500)])
def test_zoom_view_to_domain_pads_gene(hic_found, padding):
    view = gosling_spec.zoom_view_to_domain(_View(), "mm10.chr7:100000-200000", hic_found=hic_found)
    assert view.xDomain.chromosome == "chr7"
    assert view.xDomain.interval == [100000 - padding, 200000 + padding]


# ---------------------------------------------------------------- assembly arrays

class _Response:
    def __init__(self, text, status_code=200):
        self.text = text
        self.status_code = status_code


@pytest.fixture
def fake_get(monkeypatch):
    calls = []

    def install(text, status_code=200):
        def get(url, *args, **kwargs):
            calls.append(url)
            return _Response(text, status_code)
        monkeypatch.setattr(gosling_spec.requests, "get", get)
        return calls

    return install


def test_build_assembly_array(fake_get):
    calls = fake_get("chr1\t195471971\nchr2\t182113224\n")
    assert gosling_spec.build_assembly_array("mm10") == [["chr1", 195471971], ["chr2", 182113224]]
    assert calls == ["https://umgear.org/tracks/genomes/mm10/mm10.chromInfo.txt"]


def test_build_assembly_array_tolerates_blank_lines_spaces_and_three_columns(fake_get):
    fake_get("chr1\t195154279\nchrZ    82529921\nchrM\t0\t16299\n\n")
    assert gosling_spec.build_assembly_array("mm39") == [["chr1", 195154279], ["chrZ", 82529921], ["chrM", 16299]]


def test_build_assembly_array_unsupported_assembly(fake_get):
    calls = fake_get("")
    with pytest.raises(ValueError, match="not supported"):
        gosling_spec.build_assembly_array("xenLae2")
    assert calls == []


def test_build_assembly_array_http_error(fake_get):
    fake_get("Not Found", status_code=404)
    with pytest.raises(ValueError, match="Failed to retrieve"):
        gosling_spec.build_assembly_array("hg38")


def test_build_assembly_json_from_array():
    assert gosling_spec.build_assembly_json_from_array([["chr1", 100], ["chrX", 50]]) == [
        {"chrom": "chr1", "chromStart": 1, "chromEnd": 100},
        {"chrom": "chrX", "chromStart": 1, "chromEnd": 50},
    ]
    assert gosling_spec.build_assembly_json_from_array([]) == []


# ---------------------------------------------------------------- TrackSpec URL checks

@pytest.mark.parametrize("url", ["http://example.org/a.bw", "https://example.org/a.bw"])
def test_validate_url_scheme_accepts_http(url):
    gosling_spec.TrackSpec._validate_url_scheme(url)


@pytest.mark.parametrize("url", ["ftp://example.org/a.bw", "file:///a.bw", "/data/a.bw", "HTTP://example.org/a.bw"])
def test_validate_url_scheme_rejects(url):
    with pytest.raises(ValueError, match="http"):
        gosling_spec.TrackSpec._validate_url_scheme(url)


@pytest.mark.parametrize("url", ["https://example.org/a.bw", "https://example.org/A.BIGWIG", "https://example.org/a.bam.bai"])
def test_validate_url_extension_accepts(url):
    gosling_spec.TrackSpec._validate_url_extension(url, [".bw", ".bigWig", ".bam.bai"])


@pytest.mark.parametrize("url", ["https://example.org/a.bed", "https://example.org/a.bw?x=1", "https://example.org/a.bai"])
def test_validate_url_extension_rejects(url):
    with pytest.raises(ValueError, match=r"\.bw, \.bam\.bai"):
        gosling_spec.TrackSpec._validate_url_extension(url, [".bw", ".bam.bai"])
