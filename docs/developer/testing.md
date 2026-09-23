# Testing

This page describes the automated test setups that exist in the repository, how to run each one, and what CI runs. There is currently no server-side (Python/API) unit test suite; all existing tests drive the web UI in a browser.

| Setup | Location | Framework | Status |
|-------|----------|-----------|--------|
| Mocha UI tests (API mocked) | `tests/test/*.test.js` | Mocha + Playwright (`playwright`, `playwright/test` `expect`) | Active; run by CI |
| Mocha end-to-end placeholder | `tests/test/e2e/e2e.test.js` | Playwright | Stub only (comment, no tests) |
| SeleniumBase UI tests | `tests/test_*.py` | pytest + SeleniumBase | Legacy; targets v1 page markup |
| Legacy runner | `tests/run_tests`, `tests/test_commands.json`, `tests/accounts__*.py` | Selenium WebDriver | Legacy |
| Playwright scaffold | `www/js/playwright.config.js`, `www/js/tests-examples/` | `@playwright/test` | Untracked example output of `npm init playwright`; not wired up |

See also [tests/README.md](../../tests/README.md) for testing strategy notes and the list of planned UI test cases.

## Mocha + Playwright UI tests

### Layout

| File | Purpose |
|------|---------|
| `tests/package.json` | Declares `mocha` (^10.2.0), `@playwright/test` (^1.40.1), `"type": "module"`, and the script `"test": "mocha"` |
| `tests/test/helpers.js` | Shared setup: browser matrix, `gearBase` URL, `setupBrowserContext()` / `teardownBrowserContext()`, `login()` / `loginFailure()` mocks, and API mocks such as `mockGetOrganismList`, `mockGetDatasetList`, `mockGetDatasetDisplays`, `mockGetDatasetGenes`, `mockGetDatasetAnalyses`, `mockGetDatasetAvailableDisplayTypes`, `mockGetDatasetH5adInfo`, `mockGetDatasetAggregations` |
| `tests/test/gene_list_manager.test.js` | `gene_list_manager.html`: new list validation, search, facets, list actions (most complete suite) |
| `tests/test/dataset_curator.test.js` | `dataset_curator.html`: plotting, options, save (logged in / not logged in) |
| `tests/test/multigene_curator.test.js` | `multigene_curator.html`: dataset list, plotting; many `it()` bodies are still empty |
| `tests/test/compare_datasets.test.js` | Comments only (planned tests) |
| `tests/test/test_template` | Copy-and-rename template for a new page test (no `.js` extension, so Mocha ignores it) |

API responses are mocked with `page.route()` against URLs such as `${gearBase}/cgi/login.v2.cgi` and `${gearBase}/cgi/search_gene_carts.cgi`, so the tests do not need a database. The most recently registered route wins, which is how logged-in and logged-out states are switched. If a server response changes shape, update the corresponding mock in `helpers.js` or the test file.

Each suite loops over the browser list in `helpers.js` (`chromium`, `webkit`, `firefox`, `iPhone` = iPhone 12 on WebKit, `pixel` = Pixel 5 on Chromium).

### Environment variables (read in `helpers.js`)

| Variable | Effect |
|----------|--------|
| `BROWSER` | Run only the named browser (`chromium`, `webkit`, `firefox`, `iPhone`, `pixel`) |
| `LOCAL=true` | Use `http://localhost:8080` as `gearBase`; otherwise `https://devel.umgear.org` |
| `DEBUG=true` | Launch headed with `slowMo: 1000` |

### Running locally

Prerequisites: Node.js 20 (CI version) and network access to the target host (the HTML/JS is loaded from `gearBase` even though API calls are mocked).

```bash
cd tests
npm ci                                   # install mocha + @playwright/test
npx playwright install --with-deps       # download browser binaries
npm test                                 # runs mocha on ./test/*.test.js
BROWSER=chromium npm test                # single browser
BROWSER=chromium DEBUG=true npx mocha test/gene_list_manager.test.js
LOCAL=true npm test                      # against a local server on :8080
```

Mocha uses its default spec (`./test/*.{js,cjs,mjs}`, non-recursive), so `tests/test/e2e/` is not run by `npm test`. There is no `.mocharc` file; suites set `this.timeout(...)` and `this.retries(3)` themselves.

### CI

`.github/workflows/mocha_tests.yml` ("Mocha Tests") runs on push and pull request to `devel` when any `*.js` or `*.html` file changes. It uses a matrix over `chromium, firefox, webkit, iPhone, pixel` with `max-parallel: 1`, and in `./tests` it runs:

1. `npm ci`
2. `npx playwright install && npx playwright install-deps`
3. `npm run dev` in the background with `PORT=8080`
4. `npm test` with `BROWSER=<matrix value>`

The CI suite was never fully built out; filling it in is planned future work. Known gaps: `tests/package.json` defines no `dev` script, so step 3 fails silently (its output is discarded) and, because `LOCAL` is not set, the tests run against `https://devel.umgear.org`. The workflow notes a TODO to pass user/password secrets.

Other workflows in `.github/workflows/` are not test suites: `codeql.yml` (CodeQL scan of Actions, JavaScript and Python; CGIs are renamed to `.py` first), `version-js.yml` and `strip-versioning.yml` (cache-busting query strings; see [cache busting guide](./misc/cache_busting_guide.md)).

### Writing a new test

1. Copy `tests/test/test_template` to `tests/test/<page>.test.js` and replace `<PAGE>`.
2. Mock every CGI/API call the page makes with `page.route()`; add reusable mocks to `helpers.js`.
3. Prefer user-facing locators (role, text) or `data-testid` attributes over CSS selectors.
4. Do not use arrow functions for `describe` blocks that call `this.timeout()`.

## SeleniumBase / pytest tests (legacy)

Files: `tests/test_front_page.py`, `tests/test_compare_datasets.py`, `tests/test_multigene_curator.py`, `tests/test_sc_workbench.py`. Each defines `BaseCase` subclasses with `test_*` methods, reads credentials from `../gear.ini` section `[test]` (`user_email`, `password`; see [configuration](./configuration.md)), and uses SeleniumBase visual regression (`check_window`) for plots.

```bash
pip install seleniumbase pytest
cd tests                              # gear.ini is read from ../gear.ini
pytest test_front_page.py             # against https://umgear.org/
pytest test_front_page.py --data=localhost   # against http://localhost:8080/ (e.g. Docker)
```

There is no `conftest.py`, `pytest.ini` or requirements file for these tests. They target v1 selectors (`#user_email`, `#btn_sign_in`, jsTree dataset pickers) and are likely stale against the current v2 UI. `tests/.gitignore` excludes their output directories (`latest_logs`, `downloaded_files`).

## Legacy runner (`run_tests`)

`tests/run_tests` imports each module listed in `tests/test_commands.json` and calls its `main()`, which must return a list of `{"success": 0|1, "label": "..."}` dicts; it prints per-test and total pass/fail counts.

```bash
cd tests
./run_tests      # shebang is /opt/bin/python3; or: python3 run_tests
```

`test_commands.json` lists only `accounts__create_account` and `accounts__log_in`. Both use Selenium Chrome WebDriver and read `[test]` `host`, `user_name`, `user_email`, `user_institution`, `password` from `../gear.ini`. Other one-off scripts (`datasets__upload_bulk-rnaseq.py`, `gene_cart_manager__search_cart.py`, `index__primary_search__single_gene.py`) are not registered in the runner; the latter two hit `https://umgear.org/` directly and use v1 element IDs.

## Test data and environment helpers

- `tests/data/*.xlsx` - uploader spreadsheet test cases (missing sheets, count/name mismatches, non-numeric cells); described in `tests/data/test_case_descriptions.txt`.
- `tests/data/gear-test.sql` - a small MySQL dump; ignored by git (`**/*.sql` in `.gitignore`). `bin/create_test_mysql_dump.py` builds such a dump from datasets referenced by given layout IDs.
- `tests/setup_environment.py` (untracked, work in progress) - intended to create a minimal AnnData object and load `gear-test.sql` into MySQL for CI runners; the MySQL command still contains placeholder credentials.

## Playwright scaffold in `www/js/` (untracked)

`www/js/playwright.config.js`, `www/js/tests-examples/` (`demo-todo-app.spec.js`, `e2e/example.spec.js`) and `www/js/.github/workflows/playwright.yml` are the unmodified output of Playwright's project initializer. The config sets `testDir: './test/e2e'` (which does not exist under `www/js/`), five browser projects, and the HTML reporter; the nested `.github` workflow is not used by GitHub because it is not at the repository root. Treat these as reference only. If `@playwright/test`'s own runner is adopted, the config belongs in `tests/` next to `tests/test/e2e/`.

## Linting

`ruff.toml` configures Ruff for Python (includes `*.cgi`, excludes `bin/`). Run `ruff check lib www/api www/cgi listeners`. The pre-commit config (`.pre-commit-config.yaml`) currently only runs the cache-version bump hook; the Ruff hooks are not enabled.

## Related documentation

- [Code map](./code_map.md) - where page JS and CGIs live
- [API reference](./api_reference.md) - endpoints that tests mock
- [Release test plan](./misc/release_test_plan.md) - manual pre-release checklist
- [Setup guides](./setup/README.md)

---

Last updated: September 2026
