"""
Database lookups shared by the searchers that start from files on disk (projections and
cell-cell communication analyses): which datasets and gene lists the user may see, with the
details the results show.

Each returns plain dicts, so the searchers' tests can replace these functions with stubs.
"""

from gear.search.access import DATASET_LIVE_CLAUSE, dataset_visibility, gene_cart_visibility


def _placeholders(values):
    return ", ".join(["%s"] * len(values))


def access_label(user, owner_id, is_public, is_shared):
    """The access badge for a dataset: "Yours", "Public", "Shared with you" or "Group-affiliated"."""
    if user and owner_id == user.id:
        return "Yours"
    if is_public:
        return "Public"
    if is_shared:
        return "Shared with you"
    # Visible, private, and neither owned nor shared: only a group rule lets the user see it
    return "Group-affiliated"


def visible_datasets(cursor, user, ids=(), share_ids=(), ownership=(), organism_ids=()):
    """
    Datasets with one of the given ids or share_ids that the user may see, optionally narrowed by
    the ownership and organism filters. Returns a list of dicts with keys: id, share_id, title,
    ldesc, organism_id, organism, owner_id, owner_name, is_public, date_added, dtype, geo_id,
    pubmed_id, access.
    """
    ids = list(ids)
    share_ids = list(share_ids)
    if not ids and not share_ids:
        return []

    wanted = []
    params = []
    if ids:
        wanted.append(f"d.id IN ({_placeholders(ids)})")
        params.extend(ids)
    if share_ids:
        wanted.append(f"d.share_id IN ({_placeholders(share_ids)})")
        params.extend(share_ids)

    visibility, visibility_params = dataset_visibility(user, ownership)
    params.extend(visibility_params)

    conditions = ["(" + " OR ".join(wanted) + ")", DATASET_LIVE_CLAUSE, visibility]
    if organism_ids:
        conditions.append(f"d.organism_id IN ({_placeholders(organism_ids)})")
        params.extend(organism_ids)

    # Whether the dataset is shared with the user, for the access badge
    is_shared = "0"
    is_shared_params = []
    if user:
        is_shared = ("EXISTS (SELECT 1 FROM dataset_shares ds"
                     " WHERE ds.dataset_id = d.id AND ds.is_allowed = 1 AND ds.user_id = %s)")
        is_shared_params = [user.id]

    cursor.execute(
        f"""
        SELECT d.id, d.share_id, d.title, d.ldesc, d.organism_id, o.label, d.owner_id, g.user_name,
               d.is_public, d.date_added, d.dtype, d.geo_id, d.pubmed_id, {is_shared}
          FROM dataset d
               JOIN organism o ON o.id = d.organism_id
               LEFT JOIN guser g ON g.id = d.owner_id
         WHERE {" AND ".join(conditions)}
        """,
        [*is_shared_params, *params],
    )

    datasets = []
    for (dataset_id, share_id, title, ldesc, organism_id, organism, owner_id, owner_name,
         is_public, date_added, dtype, geo_id, pubmed_id, shared) in cursor.fetchall():
        datasets.append({
            "id": dataset_id,
            "share_id": share_id,
            "title": title,
            "ldesc": ldesc,
            "organism_id": organism_id,
            "organism": organism,
            "owner_id": owner_id,
            "owner_name": owner_name,
            "is_public": is_public,
            "date_added": date_added.strftime("%Y-%m-%d") if date_added else None,
            "dtype": dtype,
            "geo_id": geo_id,
            "pubmed_id": pubmed_id,
            "access": access_label(user, owner_id, is_public, shared),
        })
    return datasets


def visible_gene_carts(cursor, user, share_ids):
    """Gene lists among share_ids that the user may see: {share_id: {share_id, label, ldesc}}."""
    share_ids = list(share_ids)
    if not share_ids:
        return {}

    visibility, params = gene_cart_visibility(user, [])
    cursor.execute(
        f"SELECT gc.share_id, gc.label, gc.ldesc FROM gene_cart gc"
        f" WHERE gc.share_id IN ({_placeholders(share_ids)}) AND {visibility}",
        [*share_ids, *params],
    )
    return {share_id: {"share_id": share_id, "label": label, "ldesc": ldesc}
            for share_id, label, ldesc in cursor.fetchall()}
