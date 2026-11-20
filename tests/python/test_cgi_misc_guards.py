"""Path guards in get_dataset_display_image.cgi and the contact form checks in create_github_issue.cgi."""

import importlib.machinery
import importlib.util
import shutil
import uuid
from pathlib import Path

import pytest

from helpers.cgi_harness import CGI_DIR, run_cgi

FAKES_DIR = Path(__file__).resolve().parent / "fakes"
PREVIEWS_DIR = CGI_DIR.parent / "img" / "dataset_previews"
MISSING = "../img/dataset_previews/missing.png"


class TestGetDatasetDisplayImage:
    def run(self, dataset_id, display_id, db=None):
        return run_cgi("get_dataset_display_image.cgi", form={"dataset_id": dataset_id, "display_id": display_id},
                       db=db)

    @pytest.fixture
    def preview(self):
        dataset_id = f"pytest-{uuid.uuid4().hex}"
        path = PREVIEWS_DIR / f"{dataset_id}.5.png"
        path.write_bytes(b"png")
        yield dataset_id
        path.unlink(missing_ok=True)

    def test_existing_preview(self, preview):
        assert self.run(preview, "5").json() == f"../img/dataset_previews/{preview}.5.png"

    def test_missing_preview(self):
        assert self.run("pytest-nosuchdataset", "5").json() == MISSING

    def test_gosling_display(self, preview):
        result = self.run(preview, "5", db={"displays": [{"id": 5, "plot_type": "gosling"}]})
        assert result.json() == "../img/dataset_previews/gosling.png"

    @pytest.mark.parametrize("dataset_id, display_id", [
        ("../../cgi/get_dataset_display_image", "cgi"),
        ("/etc/passwd", "5"),
        ("x", "5/../../../../index"),
    ])
    def test_traversal_returns_missing_image(self, dataset_id, display_id):
        result = self.run(dataset_id, display_id)
        assert result.json() == MISSING
        assert "Traceback" not in result.stderr

    def test_sibling_directory_with_same_prefix_is_refused(self):
        sibling = PREVIEWS_DIR.parent / f"dataset_previews_pytest{uuid.uuid4().hex[:8]}"
        sibling.mkdir()
        try:
            (sibling / "a.5.png").write_bytes(b"png")
            result = self.run(f"../{sibling.name}/a", "5")
            assert result.json() == MISSING
        finally:
            shutil.rmtree(sibling)


def _load_github_issue():
    loader = importlib.machinery.SourceFileLoader("create_github_issue", str(CGI_DIR / "create_github_issue.cgi"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def github_issue():
    return _load_github_issue()


VALID = {"submitter_fullname": "Ann", "submitter_email": "ann@example.org", "comment_title": "Hi",
         "comment": "A question", "comment_tag": ""}


class TestValidateForm:
    def test_valid(self, github_issue):
        assert github_issue.validate_form(dict(VALID), "28") == []

    def test_required_fields(self, github_issue):
        fields = {name: "" for name in VALID}
        errors = github_issue.validate_form(fields, "28")
        assert errors == ["Name is required", "Email address is required", "Title is required",
                          "Question / Comment is required"]

    @pytest.mark.parametrize("email", ["ann", "ann@example", "a b@example.org", "ann@@example.org"])
    def test_invalid_email(self, github_issue, email):
        assert github_issue.validate_form({**VALID, "submitter_email": email}, "28") == ["Email address is not valid"]

    @pytest.mark.parametrize("name, label, limit", [
        ("comment_title", "Title", 256),
        ("comment", "Question / Comment", 60000),
        ("comment_tag", "Tags", 1000),
    ])
    def test_too_long(self, github_issue, name, label, limit):
        assert github_issue.validate_form({**VALID, name: "x" * limit}, "28") == []
        errors = github_issue.validate_form({**VALID, name: "x" * (limit + 1)}, "28")
        assert errors == [f"{label} is longer than {limit} characters"]

    @pytest.mark.parametrize("answer", ["", "27", " 28", "twenty-eight"])
    def test_wrong_security_answer(self, github_issue, answer):
        assert github_issue.validate_form(dict(VALID), answer) == ["Security check answer is incorrect"]


class TestLinkScreenshot:
    URL = "http://localhost/contact_screenshots"

    @pytest.fixture
    def root(self, github_issue, tmp_path, monkeypatch):
        monkeypatch.setattr(github_issue, "SCREENSHOT_ROOT", tmp_path)
        (tmp_path / "files").mkdir()
        (tmp_path / "files" / "shot.png").write_bytes(b"png")
        (tmp_path / "secret.png").write_bytes(b"secret")
        return tmp_path

    def links(self, root):
        return [p for p in root.iterdir() if p.is_symlink()]

    def test_uploaded_file_is_linked(self, github_issue, root):
        url = github_issue.link_screenshot("shot.png", self.URL)
        (link,) = self.links(root)
        assert url == f"{self.URL}/{link.name}" and link.suffix == ".png"
        assert str(link.readlink()) == "files/shot.png"
        assert link.read_bytes() == b"png"

    @pytest.mark.parametrize("screenshot", ["../secret.png", "nosuch.png", "", "..", "files/../../secret.png"])
    def test_names_outside_files_are_refused(self, github_issue, root, screenshot):
        url = github_issue.link_screenshot(screenshot, self.URL)
        assert url == "None (uploaded screenshot not found)"
        assert self.links(root) == []

    def test_absolute_path_is_refused(self, github_issue, root):
        assert github_issue.link_screenshot(str(root / "secret.png"), self.URL).startswith("None")
        assert self.links(root) == []

    def test_symlink_out_of_files_is_refused(self, github_issue, root):
        (root / "files" / "escape.png").symlink_to(root / "secret.png")
        assert github_issue.link_screenshot("escape.png", self.URL).startswith("None")
        assert self.links(root) == []


class TestCreateGithubIssue:
    FAKE_REQUESTS = {"requests": FAKES_DIR / "guards_requests.py",
                     "requests.exceptions": FAKES_DIR / "guards_requests.py"}

    def test_invalid_form_is_not_sent_to_github(self):
        result = run_cgi("create_github_issue.cgi", fake_modules=self.FAKE_REQUESTS,
                         form={**VALID, "submitter_email": "nope", "super_impressive_security_check": "7"})
        assert result.json() == {"error": "Email address is not valid; Security check answer is incorrect",
                                 "success": 0}
        assert not result.logged("requests.post")
