"""
SQL visibility clauses for the unified search (#1383).

Each builder returns (clause, params): a parenthesised SQL condition with %s placeholders and
the values for them. The rules match the existing search pages:

- Datasets follow www/cgi/search_datasets.cgi.
- Gene lists follow www/cgi/search_gene_carts.cgi.
- Collections (layouts) follow www/cgi/get_user_layouts.cgi.

The ownership filter is a list of "yours", "shared", "group" and "public". Several values
combine with OR. A filter with no recognized value matches nothing.
"""

# Datasets that were deleted, or never finished loading, are never searchable.
DATASET_LIVE_CLAUSE = "(d.marked_for_removal = 0 AND d.load_status = 'completed')"

# Users who share at least one group with the user given as the parameter
_GROUP_MEMBERS = (
    "SELECT DISTINCT ugm.user_id FROM user_group_membership ugm"
    " WHERE ugm.group_id IN (SELECT my.group_id FROM user_group_membership my WHERE my.user_id = %s)"
)

# Each rule is a SQL condition and whether it needs the user's id as its one parameter.
#  The order of each dict is the order the conditions appear in the clause.
_DATASET_RULES = {
    "yours": ("d.owner_id = %s", True),
    "public": ("d.is_public = 1", False),
    "shared": ("d.id IN (SELECT ds.dataset_id FROM dataset_shares ds"
               " WHERE ds.is_allowed = 1 AND ds.user_id = %s)", True),
    "group": (f"d.owner_id IN ({_GROUP_MEMBERS})", True),
}

# Sharing gene lists with a user isn't implemented, so there is no "shared" rule
_GENE_CART_RULES = {
    "yours": ("gc.user_id = %s", True),
    "public": ("gc.is_public = 1", False),
    "group": (f"gc.user_id IN ({_GROUP_MEMBERS})", True),
}

# Collections are shared by link only, so there is no "shared" rule. Site-highlighted
#  ("domain") collections are shown to everyone, so they count as public.
_LAYOUT_RULES = {
    "yours": ("l.user_id = %s", True),
    "public": ("l.is_public = 1 OR l.is_domain = 1", False),
    "group": ("l.id IN (SELECT lgm.layout_id FROM layout_group_membership lgm"
              " JOIN user_group_membership ugm ON ugm.group_id = lgm.group_id"
              " WHERE ugm.user_id = %s)", True),
}


def _build(rules, wanted, user):
    """OR together the rules named in wanted (in the rules' own order)."""
    conditions = []
    params = []
    for name, (condition, needs_user) in rules.items():
        if name not in wanted:
            continue
        conditions.append(condition)
        if needs_user:
            params.append(user.id)

    if not conditions:
        # An empty "()" would be invalid SQL; match nothing instead
        return "(0 = 1)", []
    return "(" + " OR ".join(conditions) + ")", params


def dataset_visibility(user, ownership):
    """
    Return (clause, params) limiting dataset alias "d" to what the user may see.

    Logged out: public only, and the ownership filter is ignored.
    Logged in, no filter: public, owned, or shared with the user.
    """
    if not user:
        return _build(_DATASET_RULES, {"public"}, user)
    return _build(_DATASET_RULES, set(ownership or ("public", "yours", "shared")), user)


def gene_cart_visibility(user, ownership):
    """
    Return (clause, params) limiting gene list alias "gc" to what the user may see.

    Logged out: public only. Logged in, no filter: public, owned, or owned by a group member.
    """
    if not user:
        return _build(_GENE_CART_RULES, {"public"}, user)
    return _build(_GENE_CART_RULES, set(ownership or ("public", "yours", "group")), user)


def layout_visibility(user, ownership):
    """
    Return (clause, params) limiting collection alias "l" to what the user may see.

    Logged out: public and domain only. Logged in, no filter: also owned, or in one of the
    user's groups.
    """
    if not user:
        return _build(_LAYOUT_RULES, {"public"}, user)
    return _build(_LAYOUT_RULES, set(ownership or ("public", "yours", "group")), user)
