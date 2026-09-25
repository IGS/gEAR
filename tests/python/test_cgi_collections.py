"""Dataset collection (layout) CGIs: ownership checks and single, well-formed JSON responses."""

import pytest

from helpers.cgi_harness import run_cgi

# Session "owner" is user 1, who owns collection L1; "other" is a different logged-in user.
DB = {
    "sessions": {"owner": 1, "other": 2},
    "layouts": [{"id": 10, "share_id": "L1", "user_id": 1, "label": "Old name"}],
}


class TestRenameLayout:
    def test_owner_can_rename(self):
        result = run_cgi("rename_layout.cgi", db=DB,
                         query={"session_id": "owner", "layout_share_id": "L1", "layout_name": "New name"})
        assert result.json() == {"layout_label": "New name", "layout_share_id": "L1"}
        assert result.logged("layout.save", "New name")

    def test_other_user_is_refused(self):
        result = run_cgi("rename_layout.cgi", db=DB,
                         query={"session_id": "other", "layout_share_id": "L1", "layout_name": "Hijacked"})
        assert "own" in result.json()["error"]
        assert not result.logged("layout.save")

    def test_anonymous_is_refused(self):
        result = run_cgi("rename_layout.cgi", db=DB, query={"layout_share_id": "L1", "layout_name": "X"})
        assert "logged in" in result.json()["error"]
        assert not result.logged("layout.save")

    def test_unknown_collection(self):
        result = run_cgi("rename_layout.cgi", db=DB,
                         query={"session_id": "owner", "layout_share_id": "NOPE", "layout_name": "X"})
        assert result.json()["error"] == "Dataset Collection not found."


class TestUpdateLayoutVisibility:
    def run(self, db=DB, **query):
        return run_cgi("update_layout_visibility.cgi", db=db, query=query)

    def test_owner_can_make_public(self):
        result = self.run(session_id="owner", layout_share_id="L1", visibility="true")
        assert result.json() == {"error": "", "success": 1}
        (update,) = result.writes()
        assert update["query"].startswith("UPDATE layout SET is_public")
        assert update["params"] == [1, 10]

    @pytest.mark.parametrize("session, error_fragment", [
        ("other", "collections you own"),
        (None, "must be logged in"),
    ])
    def test_non_owners_are_refused(self, session, error_fragment):
        query = {"layout_share_id": "L1", "visibility": "true"}
        if session:
            query["session_id"] = session
        result = self.run(**query)
        body = result.json()
        assert body["success"] == 0 and error_fragment in body["error"]
        assert result.writes() == []

    def test_unknown_collection(self):
        result = self.run(session_id="owner", layout_share_id="NOPE", visibility="false")
        assert result.json()["error"] == "Dataset collection not found."

    def test_database_error_prints_one_json_object(self):
        result = self.run(db={**DB, "fail_sql": ["UPDATE layout"]},
                          session_id="owner", layout_share_id="L1", visibility="false")
        body = result.json()    # fails if two JSON objects were printed
        assert body["success"] == 0 and "fake database error" in body["error"]


@pytest.mark.parametrize("script, extra", [
    ("remove_dataset_from_layout.cgi", {"dataset_id": "DS1"}),
    ("remove_display_from_layout.cgi", {"display_id": "7"}),
])
class TestRemoveFromLayout:
    def test_anonymous_gets_one_error(self, script, extra):
        result = run_cgi(script, db=DB, query={"layout_share_id": "L1", **extra})
        assert "logged in" in result.json()["error"]

    def test_unknown_collection(self, script, extra):
        result = run_cgi(script, db=DB, query={"session_id": "owner", "layout_share_id": "NOPE", **extra})
        assert result.json()["error"] == "Dataset collection not found."

    def test_other_user_is_refused(self, script, extra):
        result = run_cgi(script, db=DB, query={"session_id": "other", "layout_share_id": "L1", **extra})
        assert result.json()["success"] == 0
        assert not result.logged("layout.remove_member_by_dataset_id")
        assert not result.logged("layout.remove_member_by_display_id")

    def test_owner_removes_member(self, script, extra):
        result = run_cgi(script, db=DB, query={"session_id": "owner", "layout_share_id": "L1", **extra})
        assert result.json()["success"] == 1


class TestRemoveLayout:
    def test_anonymous(self):
        result = run_cgi("remove_layout.cgi", db=DB, query={"layout_share_id": "L1"})
        assert result.json() == {"success": 0, "error": "Not able to remove layout. User must be logged in."}

    def test_unknown_layout(self):
        result = run_cgi("remove_layout.cgi", db=DB, query={"session_id": "owner", "layout_share_id": "NOPE"})
        assert result.json()["error"] == "Not able to remove layout. Layout not found."
        assert not result.logged("layout.remove")


def test_save_layout_arrangement_requires_login():
    result = run_cgi("save_layout_arrangement.cgi", db=DB, query={"layout_share_id": "L1", "layout_arrangement": "{}"})
    assert result.json() == {"success": 0, "error": "User must be logged in"}
    assert result.writes() == []


class TestUpdateShareId:
    """update_share_id.cgi: permalink characters and the per-scope length limit (layout.share_id is VARCHAR(24))."""

    DB = {**DB, "datasets": [{"id": "D1", "share_id": "S1", "owner_id": 1}]}

    def run(self, scope, share_id, new_share_id):
        return run_cgi("update_share_id.cgi", db=self.DB, query={
            "session_id": "owner", "scope": scope, "share_id": share_id, "new_share_id": new_share_id})

    def test_collection_permalink_at_the_limit_is_saved(self):
        result = self.run("layout", "L1", "x" * 24)
        assert result.json() == {"error": "", "success": 1}
        assert result.logged("layout.save_change", "x" * 24)

    @pytest.mark.parametrize("scope, share_id, limit", [("layout", "L1", 24), ("dataset", "S1", 50)])
    def test_too_long_permalink_is_explained(self, scope, share_id, limit):
        result = self.run(scope, share_id, "x" * (limit + 1))
        body = result.json()
        assert body["success"] == 0
        assert body["error"] == f"This is not a valid permalink. It can be at most {limit} characters long."
        assert not result.logged("layout.save_change") and not result.logged("dataset.save_change")

    @pytest.mark.parametrize("new_share_id", ["my collection", "Garcia)", "_hidden", "../x"])
    def test_invalid_characters_are_explained(self, new_share_id):
        body = self.run("layout", "L1", new_share_id).json()
        assert body["success"] == 0
        assert body["error"].startswith("This is not a valid permalink. Use only letters")
