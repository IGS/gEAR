"""
fake_geardb.py - A stand-in for lib/geardb.py used when running CGI scripts in tests.

The fake is configured by a JSON "spec" file whose path is in the GEAR_FAKE_DB environment
variable, and it appends a JSON line for every SQL statement, commit and object save to the
file named by GEAR_FAKE_DB_LOG, so tests can assert what was (or was not) written.

Spec keys (all optional):
    sessions        {session_id: user_id}
    users           {user_id: {email, user_name, ...}}
    datasets        [{id, share_id, owner_id, is_downloadable, is_public, dtype, title, ...,
                      archive_path, h5ad_path, owner_display: {plotly_config}, shared_with: [user_id]}]
                    DatasetCollection.get_public/get_owned_by_user/get_shared_with_user select from these
                    (is_public defaults to 1; has_h5ad=1 means "not spatial", types filters on dtype)
    layouts         [{id, share_id, user_id, label, is_public, members: [{dataset_id, ...}]}]
    gene_carts      [{id, share_id, user_id, gctype, label}]
    genes           [{gene_symbol, organism_id, ...}]   returned by GeneCollection.get_by_gene_symbol()
    metadata        {share_id: csv text}
    displays        [{id, plot_type, ...}]   returned (as a dict) by get_display_by_id()
    verification_short   value returned by get_verification_code_short_form()
    sql             [{"match": "substring", "rows": [[...], ...] | "echo"}]
                    First matching rule wins. "echo" returns [[first param]].
    fail_sql        ["substring" | {"match": "substring", "errno": int}, ...]
                    execute() raises a database error (with that errno, if given) for matching statements
    lastrowid       value of cursor.lastrowid (default 1)

Queries against user_session are answered from "sessions" unless an sql rule matches first.
Any geardb attribute the fake doesn't provide raises AttributeError naming this file, so a
missing stub is obvious rather than silently wrong.
"""

import json
import os
import types
from pathlib import Path

# fail_sql raises a real mysql.connector.Error, so CGIs that catch only database errors handle it
try:
    from mysql.connector import Error as _DatabaseError
except ImportError:
    _DatabaseError = RuntimeError


def _database_error(msg, errno):
    """Build a database error carrying an errno, as mysql.connector would raise."""
    if _DatabaseError is RuntimeError:
        return RuntimeError(msg)
    return _DatabaseError(msg=msg, errno=errno)


_SPEC = json.loads(Path(os.environ["GEAR_FAKE_DB"]).read_text()) if os.environ.get("GEAR_FAKE_DB") else {}
_LOG = os.environ.get("GEAR_FAKE_DB_LOG")

# Like the real geardb, make json.dumps() accept gEAR objects: use a class's _serialize_json()
#  if it has one, else fall back to str() (production CGIs rely on this)
def _json_default(self, obj):
    try:
        return getattr(obj.__class__, "_serialize_json", _json_default.default)(obj)
    except Exception:
        return str(obj)


_json_default.default = json.JSONEncoder().default
json.JSONEncoder.default = _json_default

servercfg = {}
analysis_base_dir = "/tmp"
domain_url = "http://localhost"
domain_label = "gEAR"


def _log(event, **fields):
    """Record an event (SQL, commit, save) for the test to inspect."""
    if _LOG:
        with open(_LOG, "a") as fh:
            fh.write(json.dumps({"event": event, **fields}, default=str) + "\n")


class _User:
    def __init__(self, user_id):
        self.id = int(user_id)
        info = _SPEC.get("users", {}).get(str(user_id), {})
        self.email = info.get("email", f"user{user_id}@example.org")
        self.user_name = info.get("user_name", f"User {user_id}")
        self.is_admin = info.get("is_admin", 0)
        self.help_id = info.get("help_id", f"help-{user_id}")


def get_user_from_session_id(session_id=None):
    user_id = _SPEC.get("sessions", {}).get(session_id or "")
    return _User(user_id) if user_id is not None else None


def get_user_id_from_session_id(session_id=None):
    user = get_user_from_session_id(session_id)
    return user.id if user else None


def get_user_by_id(user_id=None):
    return _User(user_id) if user_id is not None else None


# ---------------------------------------------------------------- layouts (collections)

class _Layout:
    def __init__(self, spec):
        self.id = spec.get("id", 1)
        self.share_id = spec.get("share_id")
        self.user_id = spec.get("user_id")
        self.label = spec.get("label", "A collection")
        self.is_public = spec.get("is_public", 0)
        self.members = [types.SimpleNamespace(**m) for m in spec.get("members", [])]

    def load(self):
        pass

    def get_members(self):
        return self.members

    def remove_member_by_dataset_id(self, dataset_id):
        _log("layout.remove_member_by_dataset_id", layout=self.share_id, dataset_id=dataset_id)

    def remove_member_by_display_id(self, display_id):
        _log("layout.remove_member_by_display_id", layout=self.share_id, display_id=display_id)

    def remove(self):
        _log("layout.remove", layout=self.share_id)

    def save(self):
        _log("layout.save", layout=self.share_id, label=self.label)

    def save_change(self, attribute, value):
        _log("layout.save_change", layout=self.share_id, attribute=attribute, value=value)


class Layout(_Layout):
    """Fake geardb.Layout constructed directly by a CGI (e.g. add_layout.cgi); save() assigns a share ID."""

    def __init__(self, user_id=None, label=None, is_current=0, members=None, **kwargs):
        super().__init__({"id": None, "user_id": user_id, "label": label})
        self.is_current = is_current

    def save(self):
        if self.share_id is None:
            self.id, self.share_id = _SPEC.get("lastrowid", 1), "NEWLAYOUT"
        super().save()


def get_layout_by_share_id(layout_share_id):
    for spec in _SPEC.get("layouts", []):
        if spec.get("share_id") == layout_share_id:
            return _Layout(spec)
    return None


# ---------------------------------------------------------------- datasets

class Dataset(types.SimpleNamespace):
    """Fake geardb.Dataset built from a spec entry (or constructed directly by a CGI)."""

    def __init__(self, **kwargs):
        spec = next((d for d in _SPEC.get("datasets", []) if d.get("id") == kwargs.get("id")), {})
        defaults = {"is_downloadable": 1, "is_public": 1, "dtype": "single-cell-rnaseq", "title": "A dataset",
                    "ldesc": None, "pubmed_id": None, "geo_id": None, "organism_id": 1, "owner_id": None,
                    "share_id": None, "schematic_image": None, "load_status": "completed"}
        merged = {**defaults, **spec, **kwargs}
        super().__init__(**merged)

    def get_source_archive_path(self):
        return getattr(self, "archive_path", None)

    def get_tarball_path(self):
        return getattr(self, "archive_path", None) or f"/nonexistent/{self.id}.tar.gz"

    def get_file_path(self, *args, **kwargs):
        return getattr(self, "h5ad_path", None) or f"/nonexistent/{self.id}.h5ad"

    def save_change(self, attribute, value):
        _log("dataset.save_change", dataset=self.id, attribute=attribute, value=value)

    def get_owner_display(self, is_multigene=False):
        display = getattr(self, "owner_display", None)
        return types.SimpleNamespace(**display) if display else None

    def _serialize_json(self):
        return dict(vars(self))


def get_dataset_by_id(d_id=None, include_shape=None, **kwargs):
    d_id = d_id or kwargs.get("dataset_id")
    for spec in _SPEC.get("datasets", []):
        if spec.get("id") == d_id:
            return Dataset(**spec)
    return None


def get_dataset_by_share_id(share_id, *args, **kwargs):
    for spec in _SPEC.get("datasets", []):
        if spec.get("share_id") == share_id:
            return Dataset(**spec)
    return None


class DatasetCollection:
    """get_by_dataset_ids() returns no datasets; tests of search CGIs assert on the logged SQL instead."""

    def __init__(self, *args, **kwargs):
        self.datasets = []

    def get_by_dataset_ids(self, ids=None, **kwargs):
        _log("dataset_collection.get_by_dataset_ids", ids=list(ids or []))
        return self.datasets

    def _select(self, keep, dtypes):
        for spec in _SPEC.get("datasets", []):
            dtype = spec.get("dtype", "single-cell-rnaseq")
            if ((dtype in dtypes) if dtypes else dtype != "spatial") and keep(spec):
                self.datasets.append(Dataset(**spec))

    def get_public(self, has_h5ad=None, types=None, **kwargs):
        self._select(lambda d: d.get("is_public", 1) == 1, types)

    def get_owned_by_user(self, has_h5ad=None, user=None, types=None, **kwargs):
        self._select(lambda d: user is not None and d.get("owner_id") == user.id, types)

    def get_shared_with_user(self, has_h5ad=None, user=None, types=None, **kwargs):
        self._select(lambda d: user is not None and user.id in d.get("shared_with", []), types)

    def _serialize_json(self):
        return {"datasets": self.datasets}


def get_metadata_by_share_id(share_id):
    return _SPEC.get("metadata", {}).get(share_id)


def get_displays_by_user_id(user_id=None, dataset_id=None, **kwargs):
    return []


def get_display_by_id(display_id):
    for spec in _SPEC.get("displays", []):
        if str(spec.get("id")) == str(display_id):
            return dict(spec)
    return {}


def get_organism_id_by_taxon_id(taxon_id):
    return _SPEC.get("organism_id", 1)


def get_verification_code_short_form(long_form):
    return _SPEC.get("verification_short", "SHORT")


# ---------------------------------------------------------------- gene carts (lists)

class _GeneCart(types.SimpleNamespace):
    def save_change(self, attribute, value):
        _log("gene_cart.save_change", gene_cart=self.share_id, attribute=attribute, value=value)


def get_gene_cart_by_share_id(share_id, *args, **kwargs):
    for spec in _SPEC.get("gene_carts", []):
        if spec.get("share_id") == share_id:
            return _GeneCart(**spec)
    return None


# ---------------------------------------------------------------- genes

class _Gene(types.SimpleNamespace):
    def load_go_terms(self):
        pass

    def load_dbxref_links(self):
        pass

    def load_aliases(self):
        pass

    def _serialize_json(self):
        return dict(vars(self))


class GeneCollection:
    """get_by_gene_symbol() returns the spec's "genes" (unfiltered) and logs the search."""

    def __init__(self, *args, **kwargs):
        self.genes = []

    def get_by_gene_symbol(self, gene_symbol=None, exact=None, organism_id=None):
        _log("gene_collection.get_by_gene_symbol", gene_symbol=gene_symbol, exact=exact)
        self.genes = [_Gene(**g) for g in _SPEC.get("genes", [])]


# ---------------------------------------------------------------- database connection

class _Cursor:
    def __init__(self):
        self.rows = []
        self.lastrowid = _SPEC.get("lastrowid", 1)
        self.rowcount = 0

    def execute(self, query, params=()):
        flat = " ".join(str(query).split())
        params = list(params or ())
        _log("sql", query=flat, params=params)
        for rule in _SPEC.get("fail_sql", []):
            pattern, errno = (rule["match"], rule.get("errno")) if isinstance(rule, dict) else (rule, None)
            if pattern in flat:
                raise _database_error(f"fake database error for: {pattern}", errno)
        self.rows = []
        for rule in _SPEC.get("sql", []):
            if rule["match"] in flat:
                self.rows = [[params[0]]] if rule["rows"] == "echo" else [list(r) for r in rule["rows"]]
                break
        else:
            if "FROM user_session" in flat and params:
                user_id = _SPEC.get("sessions", {}).get(params[0])
                self.rows = [[user_id]] if user_id is not None else []
        self.rowcount = len(self.rows)

    def __iter__(self):
        return iter([tuple(r) for r in self.rows])

    def fetchone(self):
        return tuple(self.rows[0]) if self.rows else None

    def fetchall(self):
        return [tuple(r) for r in self.rows]

    def close(self):
        pass


class Connection:
    def __init__(self):
        pass

    def get_cursor(self, use_dict=False, **kwargs):
        return _Cursor()

    def commit(self):
        _log("commit")

    def rollback(self):
        _log("rollback")

    def close(self):
        pass


# ---------------------------------------------------------------- account, gene list and layout pages
# Used by test_cgi_account_pages.py, test_cgi_gene_list_edits.py and test_cgi_layout_pages.py. Extra spec keys:
#     servercfg          merged into servercfg, e.g. {"content": {"default_layout_share_id": "L1"}}
#     gene_carts[]       also: organism_id, ldesc, is_public, is_domain, genes (symbols), recent (1 = among
#                        the owner's recent lists), group_user_ids (users who see it as a group list)
#     layouts[]          also: is_domain, is_current, group_user_ids; members[] may have display_id,
#                        is_multigene, start_col, start_row, grid_width, grid_height
#     default_displays   [{user_id, dataset_id, is_multigene, display_id}]   for get_default_display()

servercfg.update(_SPEC.get("servercfg", {}))

_GENE_CART_DEFAULTS = {"id": None, "user_id": None, "gctype": None, "label": None, "ldesc": None,
                       "organism_id": None, "share_id": None, "is_public": None, "is_domain": None}


class Gene(_Gene):
    """Fake geardb.Gene constructed by a CGI (e.g. Gene(gene_symbol="Sox2"))."""


class _AcctGeneCart(_GeneCart):
    """A gene cart loaded by ID or by a GeneCartCollection method."""

    def add_gene(self, gene):
        self.genes.append(gene)

    def save_members(self):
        _log("gene_cart.save_members", gene_cart=self.share_id, genes=[g.gene_symbol for g in self.genes])

    def save_change(self, attribute, value):
        # Like the real GeneCart.save_change, also updates the object
        super().save_change(attribute, value)
        setattr(self, attribute, value)

    def __repr__(self):
        # Like the real GeneCart (no _serialize_json), json.dumps() turns it into a JSON string
        return json.dumps(vars(self), default=str)


def _acct_gene_cart(spec):
    return _AcctGeneCart(**{**_GENE_CART_DEFAULTS, "genes": [], **spec})


def _gene_cart_get_genes(self):
    self.genes = list(getattr(self, "genes", []))
    self.num_genes = len(self.genes)


_GeneCart.get_genes = _gene_cart_get_genes


class GeneCart(types.SimpleNamespace):
    """Fake geardb.GeneCart constructed by a CGI; save() logs its fields and assigns the spec's lastrowid."""

    def __init__(self, **kwargs):
        super().__init__(**{**_GENE_CART_DEFAULTS, **kwargs})
        if not self.share_id:
            self.share_id = f"pytest{os.urandom(4).hex()}"
        self.genes = []

    def add_gene(self, gene):
        self.genes.append(gene)

    def update_from_json(self, json_obj):
        """Same as the real GeneCart.update_from_json."""
        for key in ("label", "organism_id", "gctype", "is_public", "is_domain"):
            if key in json_obj:
                setattr(self, key, json_obj[key])
        if "session_id" in json_obj and not self.user_id:
            user = get_user_from_session_id(json_obj["session_id"])
            if user is None:
                raise Exception(f"Error: No user logged in with session_id {json_obj['session_id']}")
            self.user_id = user.id
        if "genes" in json_obj:
            self.genes = [Gene(gene_symbol=g["gene_symbol"]) for g in json_obj["genes"]]

    def save(self):
        self.id = _SPEC.get("lastrowid", 1)
        fields = {k: getattr(self, k) for k in _GENE_CART_DEFAULTS}
        _log("gene_cart.save", **fields, genes=[g.gene_symbol for g in self.genes])


def get_gene_cart_by_id(gc_id):
    for spec in _SPEC.get("gene_carts", []):
        if str(spec.get("id")) == str(gc_id):
            return _acct_gene_cart(spec)
    return None


class GeneCartCollection:
    """Selects from the spec's gene_carts; each method logs its call (with include_genes)."""

    def __init__(self, include_genes=None, **kwargs):
        self.include_genes = include_genes

    def _select(self, method, keep):
        _log(f"gene_cart_collection.{method}", include_genes=self.include_genes)
        return [_acct_gene_cart(spec) for spec in _SPEC.get("gene_carts", []) if keep(spec)]

    def get_domain(self):
        return self._select("get_domain", lambda s: s.get("is_domain") == 1)

    def get_public(self):
        return self._select("get_public", lambda s: s.get("is_public") == 1)

    def get_by_user(self, user=None):
        return self._select("get_by_user", lambda s: user is not None and s.get("user_id") == user.id)

    def get_by_user_groups(self, user=None):
        return self._select("get_by_user_groups", lambda s: user is not None and user.id in s.get("group_user_ids", []))

    def get_by_user_recent(self, user=None, n=None):
        return self._select("get_by_user_recent",
                            lambda s: user is not None and s.get("user_id") == user.id and s.get("recent"))

    def get_by_share_ids(self, share_ids=None):
        return self._select("get_by_share_ids", lambda s: s.get("share_id") in (share_ids or []))


_MULTIGENE_PLOT_TYPES = ("dotplot", "heatmap", "mg_violin", "quadrant", "volcano", "mg_tsne_static",
                         "mg_umap_static", "mg_pca_static")


class _AcctLayoutMember(types.SimpleNamespace):
    """A layout member (geardb.LayoutDisplay); save() logs the grid position it would store."""

    def __repr__(self):
        # Like the real LayoutDisplay (no _serialize_json), json.dumps() turns it into a JSON string
        return json.dumps(vars(self))

    def save(self, layout=None):
        _log("layout_member.save", layout=layout.share_id, **vars(self))


class LayoutDisplay(_AcctLayoutMember):
    """Fake geardb.LayoutDisplay; get_is_multigene() looks up the plot type in the spec's displays."""

    def __init__(self, id=None, display_id=None, grid_position=None, start_col=None, grid_width=None,
                 start_row=None, grid_height=None):
        super().__init__(id=id, display_id=display_id, dataset_id=None, is_multigene=None,
                         grid_position=grid_position, start_col=start_col, grid_width=grid_width,
                         start_row=start_row, grid_height=grid_height)

    def get_is_multigene(self):
        display = get_display_by_id(self.display_id)
        self.is_multigene = display.get("plot_type") in _MULTIGENE_PLOT_TYPES if display else None


def _layout_all_members(layout):
    """The layout's spec members as LayoutDisplay-like objects (converted once, so edits stick)."""
    if "_all_members" not in vars(layout):
        layout._all_members = [_AcctLayoutMember(**vars(m)) for m in layout.members]
        layout.members = list(layout._all_members)
    return layout._all_members


def _layout_members_where(layout, multigene):
    return [m for m in _layout_all_members(layout) if bool(getattr(m, "is_multigene", False)) == multigene]


def _layout_get_members(layout, multigene):
    layout.members = _layout_members_where(layout, multigene)
    return layout.members


def _layout_add_member(self, member):
    _log("layout.add_member", layout=self.share_id, **vars(member))


_Layout.singlegene_members = property(lambda self: _layout_members_where(self, False))
_Layout.multigene_members = property(lambda self: _layout_members_where(self, True))
_Layout.get_singlegene_members = lambda self: _layout_get_members(self, False)
_Layout.get_multigene_members = lambda self: _layout_get_members(self, True)
_Layout.add_member = _layout_add_member


class LayoutCollection:
    """Selects from the spec's layouts; each method logs its call (with include_datasets)."""

    DEFAULTS = {"id": 1, "share_id": None, "user_id": None, "label": "A collection", "is_public": 0,
                "is_domain": 0, "is_current": 0}

    def __init__(self, include_datasets=None, **kwargs):
        self.include_datasets = include_datasets

    def _select(self, method, keep):
        _log(f"layout_collection.{method}", include_datasets=self.include_datasets)
        return [types.SimpleNamespace(**{**self.DEFAULTS, **{k: v for k, v in spec.items()
                                                              if k not in ("members", "group_user_ids")}})
                for spec in _SPEC.get("layouts", []) if keep(spec)]

    def get_public(self):
        return self._select("get_public", lambda s: s.get("is_public") == 1)

    def get_domains(self):
        return self._select("get_domains", lambda s: s.get("is_domain") == 1)

    def get_by_user(self, user=None):
        return self._select("get_by_user", lambda s: user is not None and s.get("user_id") == user.id)

    def get_by_users_groups(self, user=None, append=True):
        return self._select("get_by_users_groups",
                            lambda s: user is not None and user.id in s.get("group_user_ids", []))

    def get_by_share_id(self, share_id=None):
        return self._select("get_by_share_id", lambda s: bool(share_id) and s.get("share_id") == share_id)


def get_default_display(user_id, dataset_id, is_multigene=0):
    _log("get_default_display", user_id=user_id, dataset_id=dataset_id, is_multigene=is_multigene)
    for pref in _SPEC.get("default_displays", []):
        if (pref.get("user_id") == user_id and pref.get("dataset_id") == dataset_id
                and int(pref.get("is_multigene", 0)) == int(is_multigene)):
            return pref["display_id"]
    return None


def __getattr__(name):
    raise AttributeError(
        f"fake geardb has no attribute {name!r}; add a stub for it in tests/python/helpers/fake_geardb.py"
    )


# ---------------------------------------------------------------- organisms (get_organism_list.cgi)

class OrganismCollection:
    """get_all() runs the real query, so a test supplies rows with an sql rule matching "FROM organism"."""

    def __init__(self, *args, **kwargs):
        self.organisms = []

    def get_all(self):
        cursor = Connection().get_cursor()
        cursor.execute("SELECT id, label, genus, species, strain, taxon_id FROM organism ORDER BY label")
        fields = ("id", "label", "genus", "species", "strain", "taxon_id")
        self.organisms = [types.SimpleNamespace(**dict(zip(fields, row))) for row in cursor.fetchall()]
        return self.organisms
