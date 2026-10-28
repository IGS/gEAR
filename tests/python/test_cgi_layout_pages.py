"""Dataset collection (layout) and display CGIs: listing collections and members, adding displays,
choosing a collection, and fetching or deleting saved displays."""

import json

import pytest
from helpers.cgi_harness import run_cgi

# Session "owner" is user 1; "other" is user 2.
SESSIONS = {"owner": 1, "other": 2}


def assert_clean_exit(result):
    assert result.returncode == 0, result.stderr[-2000:]
    assert "Traceback" not in result.stderr, result.stderr[-2000:]


# ---------------------------------------------------------------- get_user_layouts.cgi

class TestGetUserLayouts:
    LAYOUTS = [
        {"id": 10, "share_id": "MINE", "user_id": 1, "label": "Mine", "is_current": 1},
        {"id": 11, "share_id": "MINE2", "user_id": 1, "label": "Mine too"},
        {"id": 20, "share_id": "DOM", "user_id": 9, "label": "Hearing", "is_domain": 1, "is_public": 1,
         "is_current": 1},
        {"id": 30, "share_id": "GRP", "user_id": 3, "label": "Lab", "group_user_ids": [1]},
        {"id": 40, "share_id": "PRIV", "user_id": 3, "label": "Private"},
        {"id": 50, "share_id": "PUB", "user_id": 3, "label": "Public", "is_public": 1},
    ]
    KINDS = ["user", "domain", "group", "shared", "public"]

    def run(self, **query):
        return run_cgi("get_user_layouts.cgi", db={"sessions": SESSIONS, "layouts": self.LAYOUTS}, query=query)

    def share_ids(self, body):
        return {kind: [layout["share_id"] for layout in body[f"{kind}_layouts"]] for kind in self.KINDS}

    def test_anonymous(self):
        result = self.run()
        body = result.json()
        assert self.share_ids(body) == {"user": [], "domain": ["DOM"], "group": [], "shared": [],
                                        "public": ["DOM", "PUB"]}
        assert body["selected"] == "DOM"             # the site's current domain collection
        assert not result.logged("layout_collection.get_by_user")
        assert result.writes() == []

    def test_logged_in_user(self):
        body = self.run(session_id="owner").json()
        assert self.share_ids(body) == {"user": ["MINE", "MINE2"], "domain": ["DOM"], "group": ["GRP"], "shared": [],
                                        "public": ["DOM", "PUB"]}
        assert body["selected"] == "MINE"            # the user's current collection comes before the domain one

    def test_ownership_flag_replaces_ids(self):
        body = self.run(session_id="owner").json()
        for kind in self.KINDS:
            for layout in body[f"{kind}_layouts"]:
                assert "user_id" not in layout and "id" not in layout
                assert layout["is_owner"] == (layout["share_id"] in ("MINE", "MINE2"))

    @pytest.mark.parametrize("session", [None, "owner"])
    def test_shared_private_layout(self, session):
        query = {"layout_share_id": "PRIV"}
        if session:
            query["session_id"] = session
        body = self.run(**query).json()
        assert self.share_ids(body)["shared"] == ["PRIV"]
        assert body["selected"] == "PRIV"
        assert body["shared_layouts"][0]["is_owner"] is False

    @pytest.mark.parametrize("share_id, kind", [("MINE2", "user"), ("PUB", "public"), ("DOM", "domain"),
                                                ("GRP", "group")])
    def test_shared_layout_already_listed_is_not_repeated(self, share_id, kind):
        body = self.run(session_id="owner", layout_share_id=share_id).json()
        ids = self.share_ids(body)
        assert ids["shared"] == []
        assert share_id in ids[kind]
        assert body["selected"] == share_id

    def test_unknown_shared_layout_falls_back_to_current(self):
        body = self.run(session_id="owner", layout_share_id="NOPE").json()
        assert self.share_ids(body)["shared"] == []
        assert body["selected"] == "MINE"

    def test_no_domain(self):
        result = self.run(session_id="owner", no_domain="1")
        assert self.share_ids(result.json())["domain"] == []
        assert not result.logged("layout_collection.get_domains")

    @pytest.mark.parametrize("include_members, include_datasets", [(None, True), ("1", True), ("0", False)])
    def test_include_members(self, include_members, include_datasets):
        query = {"session_id": "owner"}
        if include_members is not None:
            query["include_members"] = include_members
        result = self.run(**query)
        calls = [e for e in result.db_log if e["event"].startswith("layout_collection.")]
        assert len(calls) == 4 and all(e["include_datasets"] is include_datasets for e in calls)


# ---------------------------------------------------------------- get_users_layout_members.cgi

class TestGetUsersLayoutMembers:
    LAYOUT = {"id": 10, "share_id": "L1", "user_id": 1, "is_public": 0, "members": [
        {"id": 1, "display_id": 101, "is_multigene": 0, "grid_position": 1},
        {"id": 2, "display_id": 201, "is_multigene": 1, "grid_position": 1},
        {"id": 3, "display_id": 102, "is_multigene": 0, "grid_position": 2},
    ]}

    def run(self, layout=None, **query):
        return run_cgi("get_users_layout_members.cgi", query=query,
                       db={"sessions": SESSIONS, "layouts": [layout or self.LAYOUT]})

    def test_owner(self):
        result = self.run(session_id="owner", layout_share_id="L1")
        body = result.json()
        # Members are JSON strings (the page JSON.parse()s each one)
        assert [json.loads(m)["display_id"] for m in body["layout_members"]["single"]] == [101, 102]
        assert [json.loads(m)["display_id"] for m in body["layout_members"]["multi"]] == [201]
        assert body["is_owner"] is True and body["is_public"] is False and body["message"] is None
        assert result.writes() == []

    @pytest.mark.parametrize("session", [None, "other"])
    def test_non_owner_sees_members(self, session):
        """A collection's share ID is enough to see its members (a shared collection overrides privacy)."""
        query = {"layout_share_id": "L1"}
        if session:
            query["session_id"] = session
        body = self.run(**query).json()
        assert body["is_owner"] is False
        assert len(body["layout_members"]["single"]) == 2 and len(body["layout_members"]["multi"]) == 1

    def test_public_layout(self):
        body = self.run(layout={**self.LAYOUT, "is_public": 1}, layout_share_id="L1").json()
        assert body["is_public"] == 1

    def test_unknown_layout(self):
        result = self.run(session_id="owner", layout_share_id="NOPE")
        assert result.json() == {"layout_members": {"single": [], "multi": []},
                                 "message": "Layout for share_id NOPE not found."}


# ---------------------------------------------------------------- add_display_to_layout.cgi

def member(id, display_id, start_col, start_row=1, grid_width=4, grid_height=1, is_multigene=0):
    return {"id": id, "display_id": display_id, "start_col": start_col, "start_row": start_row,
            "grid_width": grid_width, "grid_height": grid_height, "is_multigene": is_multigene,
            "grid_position": id}


class TestAddDisplayToLayout:
    # DS1 is user 1's private dataset; DS9 is user 9's public one; DSP is user 9's private one
    DATASETS = [{"id": "DS1", "owner_id": 1, "is_public": 0}, {"id": "DS9", "owner_id": 9, "is_public": 1},
                {"id": "DSP", "owner_id": 9, "is_public": 0}]
    DISPLAYS = [{"id": 101, "user_id": 1, "dataset_id": "DS1", "plot_type": "bar"},
                {"id": 201, "user_id": 1, "dataset_id": "DS1", "plot_type": "heatmap"},
                {"id": 301, "user_id": 9, "dataset_id": "DS9", "plot_type": "bar"},       # dataset owner's
                {"id": 302, "user_id": 3, "dataset_id": "DS9", "plot_type": "bar"},       # another user's
                {"id": 303, "user_id": 2, "dataset_id": "DS9", "plot_type": "bar"},       # user 2's own
                {"id": 401, "user_id": 9, "dataset_id": "DSP", "plot_type": "bar"},       # owner's, private dataset
                {"id": 402, "user_id": 1, "dataset_id": "DSP", "plot_type": "bar"}]       # user 1's, private dataset

    def run(self, members=(), session_id="owner", display_id="101", share_id="L1", sql=()):
        db = {"sessions": SESSIONS, "displays": self.DISPLAYS, "datasets": self.DATASETS, "sql": list(sql),
              "layouts": [{"id": 10, "share_id": "L1", "user_id": 1, "members": list(members)}]}
        query = {"layout_share_id": share_id, "display_id": display_id}
        if session_id:
            query["session_id"] = session_id
        return run_cgi("add_display_to_layout.cgi", db=db, query=query)

    def added(self, result):
        (entry,) = result.logged("layout.add_member")
        return {k: entry[k] for k in ("display_id", "grid_position", "start_col", "start_row", "grid_width",
                                      "grid_height")}

    @pytest.mark.parametrize("display_id, grid_width", [("101", 4), ("201", 6)])
    def test_add_to_empty_layout(self, display_id, grid_width):
        result = self.run(display_id=display_id)
        assert_clean_exit(result)
        assert result.json() == {"success": 1, "error": ""}
        assert self.added(result) == {"display_id": display_id, "grid_position": 1, "start_col": 1, "start_row": 1,
                                      "grid_width": grid_width, "grid_height": 1}

    def test_add_after_last_member_in_row(self):
        result = self.run(members=[member(1, 5, start_col=1), member(2, 6, start_col=5)])
        assert result.json()["success"] == 1
        assert self.added(result) == {"display_id": "101", "grid_position": 3, "start_col": 9, "start_row": 1,
                                      "grid_width": 4, "grid_height": 1}

    def test_full_row_starts_below_the_first_columns_display(self):
        """When the last row is full, the new display starts below the span of that row's first display, and
        keeps the default height (it used to take the height of that display)."""
        members = [member(1, 5, start_col=1, grid_height=2), member(2, 6, start_col=5), member(3, 7, start_col=9)]
        result = self.run(members=members)
        assert result.json()["success"] == 1
        assert self.added(result) == {"display_id": "101", "grid_position": 4, "start_col": 1, "start_row": 3,
                                      "grid_width": 4, "grid_height": 1}

    @pytest.mark.parametrize("tall_col", [5, 9])
    def test_full_row_starts_below_a_taller_display_in_any_column(self, tall_col):
        """The new row starts below the tallest display of the last row, not only the first column's."""
        members = [member(i, 4 + i, start_col=col, grid_height=3 if col == tall_col else 1)
                   for i, col in enumerate((1, 5, 9), start=1)]
        result = self.run(members=members)
        assert result.json()["success"] == 1
        assert self.added(result) == {"display_id": "101", "grid_position": 4, "start_col": 1, "start_row": 4,
                                      "grid_width": 4, "grid_height": 1}

    def test_full_row_without_a_first_column_display(self):
        """With no display in column 1 of the last row, the new row still starts below that row (it used to be
        placed on the last row itself, over its displays)."""
        members = [member(1, 5, start_col=1, grid_height=2), member(2, 6, start_col=5, start_row=2, grid_width=8)]
        result = self.run(members=members)
        assert self.added(result)["start_col"] == 1 and self.added(result)["start_row"] == 3

    def test_new_row_starts_below_a_display_from_an_earlier_row(self):
        """A display from an earlier row can reach below the last row."""
        members = [member(1, 5, start_col=1, grid_height=4), member(2, 6, start_col=5, grid_width=8),
                   member(3, 7, start_col=5, start_row=2, grid_width=8)]
        result = self.run(members=members)
        assert self.added(result)["start_col"] == 1 and self.added(result)["start_row"] == 5

    def test_space_in_last_row_covered_by_an_earlier_display(self):
        """The free columns at the end of the last row are used only if no earlier tall display covers them."""
        members = [member(1, 5, start_col=1), member(2, 6, start_col=5), member(3, 7, start_col=9, grid_height=3),
                   member(4, 8, start_col=1, start_row=2, grid_width=8)]
        result = self.run(members=members)
        assert self.added(result)["start_col"] == 1 and self.added(result)["start_row"] == 4

    def test_display_that_would_overflow_the_row_starts_a_new_row(self):
        """A display that doesn't fit in the columns left in the last row goes on a new row (it used to be
        placed past column 12)."""
        result = self.run(members=[member(1, 5, start_col=1, grid_width=9)])
        assert self.added(result)["start_col"] == 1 and self.added(result)["start_row"] == 2

    def test_multigene_members_are_placed_separately(self):
        result = self.run(members=[member(1, 5, start_col=1), member(2, 6, start_col=5)], display_id="201")
        assert self.added(result) == {"display_id": "201", "grid_position": 1, "start_col": 1, "start_row": 1,
                                      "grid_width": 6, "grid_height": 1}

    def test_legacy_layout_is_rebuilt(self):
        """A layout where every member has start_col 1 (from before the grid) gets real grid positions first."""
        members = [member(1, 5, start_col=1), member(2, 6, start_col=1), member(3, 7, start_col=1),
                   member(4, 8, start_col=1)]
        result = self.run(members=members)
        assert result.json()["success"] == 1
        saves = [(e["id"], e["start_col"], e["start_row"]) for e in result.logged("layout_member.save")]
        assert saves == [(1, 1, 1), (2, 5, 1), (3, 9, 1), (4, 1, 2)]
        assert self.added(result)["start_col"] == 5 and self.added(result)["start_row"] == 2

    @pytest.mark.parametrize("session, error", [(None, "must be logged in"), ("expired", "must be logged in"),
                                                ("other", "doesn't own it")])
    def test_non_owner_is_refused(self, session, error):
        legacy = [member(1, 5, start_col=1), member(2, 6, start_col=1)]
        result = self.run(members=legacy, session_id=session)
        body = result.json()
        assert body["success"] == 0 and error in body["error"]
        assert not result.logged("layout.add_member")
        assert not result.logged("layout_member.save")      # the legacy rebuild doesn't run either
        assert result.writes() == []

    @pytest.mark.parametrize("session, display_id", [("owner", "301"), ("other", "303")])
    def test_dataset_owners_or_own_display_on_public_dataset(self, session, display_id):
        """A display saved by the dataset's owner can be added to anyone's layout, as can the user's own."""
        layout_owner = {"owner": 1, "other": 2}[session]
        db = {"sessions": SESSIONS, "displays": self.DISPLAYS, "datasets": self.DATASETS,
              "layouts": [{"id": 10, "share_id": "L1", "user_id": layout_owner}]}
        result = run_cgi("add_display_to_layout.cgi", db=db,
                         query={"session_id": session, "layout_share_id": "L1", "display_id": display_id})
        assert result.json() == {"success": 1, "error": ""}
        assert self.added(result)["display_id"] == display_id

    def test_dataset_owners_display_on_dataset_shared_with_user(self):
        result = self.run(display_id="401", sql=[{"match": "FROM dataset_shares", "rows": [[1]]}])
        assert result.json() == {"success": 1, "error": ""}
        (share_check,) = result.logged("sql", "FROM dataset_shares")
        assert share_check["params"] == ["DSP", 1]
        assert "is_allowed = 1" in share_check["query"]

    @pytest.mark.parametrize("display_id, error", [
        ("302", "belongs to another user"),            # not the user's own nor the dataset owner's
        ("401", "don't have access to its dataset"),   # the dataset owner's, but the dataset is private
        ("402", "don't have access to its dataset"),   # the user's own, but the dataset is no longer shared
        ("999", "Display not found"),
        ("", "Display not found"),
    ])
    def test_display_the_user_may_not_use_is_refused(self, display_id, error):
        legacy = [member(1, 5, start_col=1), member(2, 6, start_col=1)]
        result = self.run(members=legacy, display_id=display_id)
        assert_clean_exit(result)
        body = result.json()
        assert body["success"] == 0 and error in body["error"]
        assert not result.logged("layout.add_member")
        assert not result.logged("layout_member.save")      # the legacy rebuild doesn't run either
        assert result.writes() == []

    def test_unknown_layout(self):
        result = self.run(share_id="NOPE")
        assert result.json() == {"success": 0, "error": "Layout not found"}
        assert not result.logged("layout.add_member")


# ---------------------------------------------------------------- save_user_chosen_layout.cgi

class TestSaveUserChosenLayout:
    LOOKUP = "FROM layout WHERE share_id"
    GROUPS = "FROM layout_group_membership"

    def run(self, layout_row=(10, 1, 0, 0), in_group=False, **query):
        """layout_row: (id, user_id, is_public, is_domain) of the looked-up layout, or None if unknown."""
        sql = [{"match": self.LOOKUP, "rows": [list(layout_row)] if layout_row else []}]
        if in_group:
            sql.append({"match": self.GROUPS, "rows": [[1]]})
        return run_cgi("save_user_chosen_layout.cgi", db={"sessions": SESSIONS, "sql": sql}, query=query)

    def test_logged_in_user_saves_choice(self):
        result = self.run(session_id="owner", layout_share_id="L1")
        assert result.json() == {"success": 1, "message": ""}
        (lookup,) = result.logged("sql", self.LOOKUP)
        assert lookup["params"] == ["L1"]
        (update,) = result.writes()
        assert update["query"] == "UPDATE guser SET layout_id = %s WHERE id = %s"
        assert update["params"] == [10, 1]
        assert result.logged("commit")

    @pytest.mark.parametrize("layout_row", [(20, 9, 1, 0), (21, 9, 0, 1), (22, 9, 1, 1)])
    def test_public_or_site_layout_of_another_user(self, layout_row):
        """The collection selectors list public and site (domain) collections, so they can be chosen."""
        result = self.run(layout_row=layout_row, session_id="owner", layout_share_id="PUB")
        assert result.json() == {"success": 1, "message": ""}
        (update,) = result.writes()
        assert update["params"] == [layout_row[0], 1]

    def test_private_layout_in_users_group(self):
        """The collection selectors list collections shared with the user's groups, so they can be chosen."""
        result = self.run(layout_row=(30, 3, 0, 0), in_group=True, session_id="owner", layout_share_id="GRP")
        assert result.json() == {"success": 1, "message": ""}
        (group_check,) = result.logged("sql", self.GROUPS)
        assert group_check["params"] == [30, 1]
        (update,) = result.writes()
        assert update["params"] == [30, 1]

    def test_another_users_private_layout_is_refused(self):
        """Another user's private collection can't be made the user's chosen one, even with its share ID."""
        result = self.run(layout_row=(40, 3, 0, 0), session_id="owner", layout_share_id="PRIV")
        body = result.json()
        assert body["success"] == 0 and "belongs to another user" in body["message"]
        assert result.writes() == []
        assert not result.logged("commit")

    def test_unknown_layout_clears_the_choice(self):
        """An unknown share ID stores no layout, so get_session_info.v2.cgi falls back to the site default."""
        result = self.run(layout_row=None, session_id="owner", layout_share_id="NOPE")
        assert result.json()["success"] == 1
        (update,) = result.writes()
        assert update["params"] == [None, 1]

    @pytest.mark.parametrize("session", [None, "expired"])
    def test_anonymous_is_refused(self, session):
        query = {"layout_share_id": "L1"}
        if session:
            query["session_id"] = session
        result = self.run(**query)
        assert result.json() == {"success": 0, "message": "Invalid session_id"}
        assert result.writes() == []


# ---------------------------------------------------------------- get_default_display.cgi

class TestGetDefaultDisplay:
    # Dataset DS1 is owned by user 9, who has set displays 900 (single-gene) and 901 (multigene).
    # Display 900 has been placed in saved collections, so it has a suggested auto-layout width.
    DB = {"sessions": SESSIONS, "datasets": [{"id": "DS1", "owner_id": 9}],
          "default_displays": [{"user_id": 9, "dataset_id": "DS1", "is_multigene": 0, "display_id": 900},
                               {"user_id": 9, "dataset_id": "DS1", "is_multigene": 1, "display_id": 901},
                               {"user_id": 1, "dataset_id": "DS1", "is_multigene": 0, "display_id": 100}],
          "auto_layout_widths": [{"display_id": 900, "width": 6, "rule": "widest width used by >= 25% of uses (2/3)",
                                  "width_counts": {6: 2, 4: 1}}]}

    NO_WIDTH = {"grid_width": None, "width_rule": "display is in no saved collection", "width_counts": {}}

    def run(self, **query):
        return run_cgi("get_default_display.cgi", db=self.DB, query=query)

    def test_users_own_default(self):
        result = self.run(session_id="owner", dataset_id="DS1")
        assert result.json() == {"default_display_id": 100, **self.NO_WIDTH}
        assert [e["user_id"] for e in result.logged("get_default_display")] == [1]
        assert [e["display_id"] for e in result.logged("get_auto_layout_width")] == [100]
        assert result.writes() == []

    @pytest.mark.parametrize("session", [None, "other"])
    def test_falls_back_to_dataset_owners_default(self, session):
        query = {"dataset_id": "DS1"}
        if session:
            query["session_id"] = session
        assert self.run(**query).json() == {
            "default_display_id": 900, "grid_width": 6,
            "width_rule": "widest width used by >= 25% of uses (2/3)", "width_counts": {"6": 2, "4": 1}}

    def test_multigene(self):
        result = self.run(session_id="owner", dataset_id="DS1", is_multigene="1")
        assert result.json() == {"default_display_id": 901, **self.NO_WIDTH}

    def test_unknown_dataset(self):
        result = self.run(session_id="other", dataset_id="NOPE")
        assert_clean_exit(result)
        assert result.json() == {"default_display_id": None, "grid_width": None, "width_rule": None, "width_counts": {}}
        assert result.logged("get_auto_layout_width") == []


# ---------------------------------------------------------------- get_dataset_display.cgi

class TestGetDatasetDisplay:
    DISPLAY = {"id": 101, "dataset_id": "DS1", "user_id": 1, "label": "Bar", "plot_type": "bar",
               "plotly_config": {"gene_symbol": "Sox2"}}

    def test_known_display(self):
        result = run_cgi("get_dataset_display.cgi", db={"displays": [self.DISPLAY]}, query={"display_id": "101"})
        assert result.json() == self.DISPLAY

    @pytest.mark.parametrize("query", [{"display_id": "999"}, {}])
    def test_unknown_display(self, query):
        result = run_cgi("get_dataset_display.cgi", db={"displays": [self.DISPLAY]}, query=query)
        assert result.json() == {}


# ---------------------------------------------------------------- delete_dataset_display.cgi

class TestDeleteDatasetDisplay:
    DB = {"sessions": SESSIONS, "displays": [{"id": 101, "dataset_id": "DS1", "user_id": 1, "plot_type": "bar"}]}

    def run(self, **query):
        return run_cgi("delete_dataset_display.cgi", db=self.DB, query=query)

    def test_owner_deletes(self):
        result = self.run(session_id="owner", id="101")
        assert result.json() == {"success": True}
        (delete,) = result.writes()
        assert delete["query"] == "DELETE FROM dataset_display where id = %s"
        assert delete["params"] == ["101"]
        assert result.logged("commit")

    @pytest.mark.parametrize("query", [
        {"session_id": "other", "id": "101"},
        {"session_id": "expired", "id": "101"},
        {"id": "101"},
        {"session_id": "owner"},
    ])
    def test_refused(self, query):
        result = self.run(**query)
        assert_clean_exit(result)
        assert result.json() == {"success": False}
        assert result.writes() == []

    def test_unknown_display_counts_as_deleted(self):
        result = self.run(session_id="owner", id="999")
        assert result.json() == {"success": True}
        assert result.writes() == []
