# Testing

This page covers the automated tests in the repository: the server-side pytest suite, which runs in CI, and the older UI suites kept for reference.

| Suite | Location | Runs in CI | Needs a live site or MySQL? |
|-------|----------|------------|------------------------------|
| **Server-side pytest suite** | `tests/python/` | Yes: `python_tests.yml` | No |
| Mocha + Playwright UI tests | `tests/legacy/test/` | Manual only (`mocha_tests.yml`) | Yes (loads pages from `devel.umgear.org` or `localhost:8080`) |
| SeleniumBase UI tests (v1 UI) | `tests/legacy/test_*.py` | No | Yes |
| Selenium runner (v1 UI) | `tests/legacy/run_tests` | No | Yes |

## Server-side pytest suite (`tests/python`)

The suite covers:

- **CGI scripts**, run against a fake database.
- **Flask API resources**, through Flask's test client.
- **Upload processors** (`gear.anndata_processor`, `gear.spatial_processor`, `gear.spatialhandler`), using the real example files in `www/user_templates/`.
- **Guards:** every server-side Python file parses, `openapi.yaml` validates, and relative Markdown links resolve.

### Layout

| Path | Purpose |
|------|---------|
| `conftest.py` | Puts `lib/` and `www/api/` on `sys.path`. Replaces `gear.db.MySQLDB.connect` with a fake before `geardb` is imported (see note below). Fixtures: `upload_session`, `tmp_cart`, `example_mex_tar`, `example_3tab_tar_gz`. |
| `helpers/cgi_harness.py` | `run_cgi(...)` runs a CGI in a subprocess and returns a `CGIResult` |
| `helpers/run_cgi.py` | Runner used by `run_cgi`: installs the fake `geardb`, then runs the script as `__main__` |
| `helpers/fake_geardb.py` | Spec-driven stand-in for `lib/geardb.py` that logs every SQL statement, commit and save |
| `fakes/` | Small fake modules for individual tests, such as a `gear.analysis` that returns a synthetic AnnData |
| `test_cgi_*.py` | CGI tests, grouped by area (uploads, collections, datasets, accounts, downloads, gene lists and search, comparison) |
| `test_api_*.py` | Flask resources (`Aggregations`, `TrackHubCopy`) |
| `test_anndata_processor.py`, `test_spatial.py` | Upload processing, including cancellation when an upload is deleted mid-job |
| `test_geardb_helpers.py` | Helpers in the real `geardb` module |
| `test_syntax.py`, `test_docs.py` | Guards |
| `requirements.txt`, `requirements-spatial.txt` | Pinned to the production versions in `docker/requirements.txt` |

> **Note:** importing `geardb` opens a MySQL connection, because `LayoutCollection._cnx = Connection()` runs when the class is defined. `conftest.py` patches the connection so the real module can be imported without a database or `gear.ini`.

### Running locally

Use Python 3.14, the production version:

```bash
python3.14 -m venv ~/venvs/gear-tests && source ~/venvs/gear-tests/bin/activate
pip install -r tests/python/requirements.txt          # or requirements-spatial.txt for the spatial tests
cd tests/python
python -m pytest                    # all tests; spatial ones are skipped if spatialdata isn't installed
python -m pytest -m "not spatial"   # the "core" CI job
python -m pytest test_cgi_downloads.py -k owner -v
```

The CGI tests create and remove their own directories under `www/uploads/files/` and files under `www/carts/`; both folders are gitignored.

### How CGI tests work

`run_cgi()` runs `www/cgi/<script>` in a subprocess from `www/cgi` (many scripts find `lib/` relative to that directory):

- **Parameters:** it sets `REQUEST_METHOD` and `QUERY_STRING` for `query=`, or a multipart body for `form=`/`files=`, plus `HTTP_COOKIE` for `cookies=`.
- **Fake database:** the runner puts `helpers/fake_geardb.py` into `sys.modules` as `geardb` before the script runs. That way it's used even by scripts that put the real `lib/` first on `sys.path`.
- **Other fakes:** other modules can be replaced with `fake_modules={"gear.analysis": path}`.

The fake database is configured with a JSON-style spec:

```python
from helpers.cgi_harness import run_cgi

DB = {
    "sessions": {"owner": 1, "other": 2},                       # session_id -> user id
    "layouts": [{"id": 10, "share_id": "L1", "user_id": 1}],
    "datasets": [{"id": "DS1", "share_id": "S1", "owner_id": 1, "is_downloadable": 0}],
    "sql": [{"match": "FROM note", "rows": [[1]]}],            # canned rows for matching queries
    "fail_sql": ["UPDATE layout"],                              # make matching statements raise
}

def test_other_user_cannot_rename():
    result = run_cgi("rename_layout.cgi", db=DB,
                     query={"session_id": "other", "layout_share_id": "L1", "layout_name": "X"})
    assert "own" in result.json()["error"]     # .json() also asserts exactly one JSON document
    assert not result.logged("layout.save")    # nothing was saved
    assert result.writes() == []               # no INSERT/UPDATE/DELETE was executed
```

`CGIResult` also has:

- `.status`: the HTTP status from a `Status:` header, 200 by default
- `.headers` (lower-cased names), `.body` and `.text`
- `.json_documents()`: every JSON document printed

If a script needs a `geardb` function the fake doesn't have, the fake raises an `AttributeError` naming `fake_geardb.py`. Add a small stub there.

### How in-process tests work

Tests of library code and Flask resources import the real modules directly, because `conftest.py` has already patched the database connection. They use pytest's `monkeypatch` to replace what the code under test would read from the database or disk:

- `test_anndata_processor.py` points `gear.anndata_processor.UPLOADS_BASE_DIR` at `tmp_path` and processes copies of the example archives.
- `test_api_aggregations.py` replaces `get_adata_*` with a small AnnData and calls the resource through `Flask.test_client()`.

Import individual resources (`from resources.aggregations import Aggregations`), not `www/api/api.py`: that module calls `setrlimit` and imports `gear.orthology`, which queries the database at import.

### Markers

`spatial` marks tests that need the SpatialData stack (`requirements-spatial.txt`). They're skipped automatically when `spatialdata` isn't installed.

### CI

`.github/workflows/python_tests.yml` runs on push and pull request to `devel` and `main` when server-side code, docs, example files or the tests change. It can also be started by hand (`workflow_dispatch`). It has two jobs on Python 3.14:

- **core:** `pip install -r tests/python/requirements.txt`, then `pytest -m "not spatial"`.
- **spatial:** installs `requirements-spatial.txt`, then runs `pytest -m spatial`.

Other workflows are not test suites:

- `codeql.yml`: CodeQL scan of Actions, JavaScript and Python. CGIs are renamed to `.py` first.
- `version-js.yml` and `strip-versioning.yml`: cache-busting query strings; see the [cache busting guide](./misc/cache_busting_guide.md).

### Adding tests

- **Changing a CGI or API response:** add or update a test that covers both the success path and the error you're guarding against. Assert on `.json()`, so a script that prints two JSON objects or an empty body fails.
- **Permission checks:** always assert that a refused request made no writes (`result.writes() == []`).
- **Processors:** build small inputs in `tmp_path` or reuse the files in `www/user_templates/`. Keep each test to a few seconds.

## Legacy UI suites (`tests/legacy`)

These suites are not maintained. They're kept for reference, especially the notes on mocking with Playwright's `page.route()` and the list of planned UI tests in [tests/legacy/README.md](../../tests/legacy/README.md).

- **Mocha + Playwright** (`tests/legacy/test/*.test.js`): run with `cd tests/legacy && npm ci && npx playwright install && npm test`. API calls are mocked, but the pages are loaded from `https://devel.umgear.org`, or `http://localhost:8080` with `LOCAL=true`. `mocha_tests.yml` only runs when started by hand.
- **SeleniumBase** (`tests/legacy/test_*.py`) and the **Selenium runner** (`tests/legacy/run_tests`): written for the v1 UI and run against a live site, reading credentials from `gear.ini` `[test]`.

`www/js/playwright.config.js` and `www/js/tests-examples/` (untracked) are the unmodified output of `npm init playwright` and are not wired up.

## Test data

- `tests/data/*.xlsx`: uploader spreadsheet test cases, described in `tests/data/test_case_descriptions.txt`.
- `tests/data/gear-test.sql`: a small MySQL dump, gitignored. `bin/create_test_mysql_dump.py` builds one from given layouts.
- `www/user_templates/example_mex.tar` and `example_3tab.tar.gz`: example uploads, also used by the pytest suite.

## Linting

`ruff.toml` configures Ruff for Python (it includes `*.cgi` and excludes `bin/`). Run `ruff check lib www/api www/cgi listeners`. Ruff isn't run in CI, and the pre-commit config only runs the cache-version bump hook.

## Related documentation

- [Code map](./code_map.md): where page JS and CGIs live
- [API reference](./api_reference.md) and [OpenAPI spec](./openapi.yaml): the endpoints under test
- [Release test plan](./misc/release_test_plan.md): manual pre-release checklist
- [Setup guides](./setup/README.md)

---

Last updated: September 2026
