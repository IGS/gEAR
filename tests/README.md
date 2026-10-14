# Testing

| Directory | What it is |
| --- | --- |
| [`python/`](python) | **Maintained, run in CI.** Server-side pytest suite: CGI scripts (against a fake database), Flask API resources, the upload processors, and documentation checks. No live server, MySQL or `gear.ini` needed. |
| [`legacy/`](legacy) | Older UI suites (Mocha + Playwright, SeleniumBase, a Selenium runner) that run against a live site. Kept for reference; see [`legacy/README.md`](legacy/README.md). |
| [`data/`](data) | Fixture files (metadata spreadsheets, `gear-test.sql`). |

## Running the Python suite

Use Python 3.14 (the production version):

```bash
python3.14 -m venv ~/venvs/gear-tests
source ~/venvs/gear-tests/bin/activate
pip install -r tests/python/requirements.txt           # add requirements-spatial.txt for spatial tests

cd tests/python
python -m pytest                     # everything; what the "core" CI job runs (spatial tests skip without spatialdata)
python -m pytest -m spatial          # only tests marked "spatial" (markers are explained in pytest.ini)
python -m pytest test_cgi_uploads.py -k tar
```

CI runs this suite in [`.github/workflows/python_tests.yml`](../.github/workflows/python_tests.yml) on pushes and pull requests to `devel` and `main`.

The [testing guide](../docs/developer/testing.md) explains how the fake database and CGI harness work and how to add tests.
