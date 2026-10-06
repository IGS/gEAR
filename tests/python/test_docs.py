"""Documentation checks: the OpenAPI spec validates and relative Markdown links resolve."""

import re
from pathlib import Path

import pytest
import yaml
from openapi_spec_validator import validate

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS = REPO_ROOT / "docs"
# Historical material that is kept as-is
SKIP_DIRS = {"posters", "ui-v2-design"}


def test_openapi_spec_is_valid():
    validate(yaml.safe_load((DOCS / "developer" / "openapi.yaml").read_text()))


def _slug(heading: str) -> str:
    """GitHub-style heading anchor."""
    heading = re.sub(r"<[^>]+>", "", heading).strip().lower()
    return re.sub(r"[^\w\- ]", "", heading).replace(" ", "-")


def _anchors(md_file: Path) -> set:
    anchors, seen = set(), {}
    for line in md_file.read_text(errors="ignore").splitlines():
        match = re.match(r"^#{1,6}\s+(.*?)\s*#*$", line)
        if match:
            slug = _slug(match.group(1))
            count = seen.get(slug, 0)
            seen[slug] = count + 1
            anchors.add(slug if count == 0 else f"{slug}-{count}")
        anchors.update(re.findall(r'<a (?:name|id)="([^"]+)"', line))
    return anchors


def _markdown_files():
    files = [p for p in DOCS.rglob("*.md") if not SKIP_DIRS & set(p.parts)]
    return sorted(files + [REPO_ROOT / "README.md", REPO_ROOT / "tests" / "README.md"])


def _links(md_file: Path):
    text = md_file.read_text(errors="ignore")
    text = re.sub(r"```.*?```", "", text, flags=re.S)   # fenced code
    text = re.sub(r"`[^`\n]*`", "", text)                # inline code
    yield from re.findall(r"\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)", text)
    yield from re.findall(r'src="([^"]+)"', text)


@pytest.mark.parametrize("md_file", _markdown_files(), ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_relative_links_resolve(md_file):
    broken = []
    for link in _links(md_file):
        if re.match(r"^(https?:|mailto:|data:)", link):
            continue
        path, _, anchor = link.partition("#")
        target = md_file if not path else (md_file.parent / path).resolve()
        if not target.exists():
            broken.append(f"missing file: {link}")
        elif anchor and target.suffix == ".md" and anchor.lower() not in _anchors(target):
            broken.append(f"missing anchor: {link}")
    assert not broken, "\n".join(broken)
