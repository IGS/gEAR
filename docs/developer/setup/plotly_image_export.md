# Plotly image export (kaleido and Chromium)

gEAR renders Plotly figures to static images on the server in a few places:

- the "Download image" menu item on Plotly tiles (`www/api/resources/plotly_data.py`, `www/api/resources/mg_plotly_data.py`)
- the preview thumbnail made when a Plotly display is saved (`www/cgi/save_dataset_display.cgi`)
- batch preview thumbnails (`bin/generate_static_display_images.py`)
- the image downloads on the cell-cell communication page (`www/cgi/get_domino_plots.cgi`)

Plotly does this with the [kaleido](https://github.com/plotly/Kaleido) package. Since version 1.0, kaleido does not include a browser: it drives a separately installed Chrome or Chromium. So every web server needs three things:

1. the `kaleido` Python package
2. a Chrome or Chromium browser that the Apache user (`www-data`) can run, plus the system libraries it needs
3. the `sys.executable` fix in `www/api/api.wsgi` (already in the repository; see [Why api.wsgi sets sys.executable](#why-apiwsgi-sets-sysexecutable))

When any of these is missing, "Download image" on a Plotly tile fails with an HTTP 500 and the user sees "The image could not be downloaded". The cause is in the Apache error log; see [Troubleshooting](#troubleshooting).

## Server setup

These steps are written for Ubuntu (the Chromium install was tested on 22.04). They assume Python is installed at `/opt/Python-${PYTHONV}` as described in [the Python setup](./python.md).

```bash
export PYTHONV=3.14.4
```

### 1. Install kaleido

`kaleido` is in `docker/requirements.txt`, so a server set up from that file already has it. To check, or to install it on an existing server:

```bash
/opt/Python-${PYTHONV}/bin/python3 -c "import kaleido; print(kaleido.__version__)"
sudo /opt/Python-${PYTHONV}/bin/python3 -m pip install kaleido==1.4.0
```

Keep the version in step with `docker/requirements.txt`.

### 2. Install Chrome in a shared location

plotly includes a command, `plotly_get_chrome`, that downloads a build of Chrome (Chrome for Testing) that kaleido is tested with. By default it installs into the home directory of the user who runs it (with `sudo`, under `/root/.local/share/`), which the Apache user cannot see. So give it a shared directory with `--path`, and link the browser into `/usr/bin`:

```bash
sudo mkdir -p /opt/chrome
sudo /opt/Python-${PYTHONV}/bin/plotly_get_chrome -y --path /opt/chrome   # prints the Chrome path
sudo ln -sf /opt/chrome/chrome-linux64/chrome /usr/bin/google-chrome
sudo chmod -R a+rX /opt/chrome
```

If the printed path differs from `/opt/chrome/chrome-linux64/chrome`, use the printed one in the `ln` command. `/usr/bin/google-chrome` is one of the paths kaleido searches, so the link works no matter which user or `PATH` Apache runs with.

Chrome also needs system libraries that a server usually does not have. List the missing ones:

```bash
ldd /opt/chrome/chrome-linux64/chrome | grep "not found"
```

Then install the packages that provide them. On Ubuntu 22.04 these were:

```bash
sudo apt-get update
sudo apt-get install libatk1.0-0 libatk-bridge2.0-0 libcups2 libxkbcommon0 libatspi2.0-0 libgbm1 libasound2
```

On Ubuntu 24.04, several packages have a `t64` suffix:

```bash
sudo apt-get install libatk1.0-0t64 libatk-bridge2.0-0t64 libcups2t64 libxkbcommon0 libatspi2.0-0t64 libgbm1 libasound2t64
```

Run the `ldd` command again until it prints nothing. If it lists a library not covered here, `apt-file search <library>.so` finds the package. Then check that Chrome starts as the Apache user. It should print `<html><head></head><body><p>ok</p></body></html>`:

```bash
sudo -u www-data HOME=/var/www /usr/bin/google-chrome --headless --no-sandbox --disable-gpu \
  --dump-dom 'data:text/html,<p>ok</p>'
```

Leave out `--headless` and Chrome tries to open a window and exits with "Missing X server or $DISPLAY", which says nothing about whether it works.

Chrome for Testing is published for amd64 only. For an arm64 server, or to have the system libraries installed for you, use [Playwright's Chromium](#alternative-playwright-chromium) instead.

### 3. Check which browser kaleido uses

kaleido uses the first browser it finds, in this order:

1. a browser previously downloaded by `plotly_get_chrome` (or `kaleido_get_chrome`, `choreo_get_chrome`) into the Apache user's own home directory
2. Google Chrome: on the `PATH`, then `/usr/bin/google-chrome-stable`, `/usr/bin/google-chrome` and `/usr/bin/chrome`
3. Chromium: on the `PATH`, then `/usr/bin/chromium` and `/usr/bin/chromium-browser`

So a leftover browser earlier in that order wins over the one you just installed, even if it does not work. List what is there:

```bash
ls -l /usr/bin/google-chrome* /usr/bin/chrome /usr/bin/chromium* 2>/dev/null
ls -d /var/www/.local/share/choreographer 2>/dev/null
```

Remove anything that is not the browser from step 2.

### 4. Update the code and restart Apache

Pull a version of gEAR that includes the `sys.executable` fix in `www/api/api.wsgi`, then restart Apache so the API processes pick up kaleido and the fix:

```bash
sudo systemctl restart apache2
```

### 5. Check it

First, as the Apache user from a shell. This prints the size of a small PDF:

```bash
sudo -u www-data HOME=/var/www /opt/Python-${PYTHONV}/bin/python3 -c \
  "import plotly.graph_objects as go; print(len(go.Figure().to_image(format='pdf')), 'bytes')"
```

Then in the browser: open a dataset with a Plotly display (bar, violin, scatter, etc.) and use "Download image" in the tile's menu. It downloads the metadata CSV and a PDF.

Do both checks. The shell check passing does not prove the download works: Apache runs the API differently from a shell (see the next section).

## Why api.wsgi sets sys.executable

kaleido's browser library (choreographer) does not start Chromium directly. It runs `sys.executable` with a small Python helper script, and the helper starts Chromium. The helper imports packages that kaleido installed into `/opt/Python-${PYTHONV}`.

Under mod_wsgi, `sys.executable` is not always the interpreter running the API. On one production server it was the system `/usr/bin/python3`. The helper failed on `import logistro` and exited before Chromium started. The Apache error log only said "The browser seemed to close immediately after starting", and every test run from a shell passed.

`www/api/api.wsgi` now points `sys.executable` at `sys.prefix/bin/pythonX.Y` before loading the app. With `WSGIPythonHome "/opt/Python-${PYTHONV}"` (see [the Apache setup](./apache.md)), that is `/opt/Python-${PYTHONV}/bin/python3.14`.

To see what mod_wsgi reports on a server, temporarily add this line after `import sys` in `api.wsgi`, restart Apache and load any page:

```python
print(f"api.wsgi sys.prefix={sys.prefix!r} sys.executable={sys.executable!r}", file=sys.stderr)
```

```bash
sudo grep "api.wsgi sys.prefix" /var/log/apache2/*error.log | tail -3
```

CGI scripts such as `save_dataset_display.cgi` are not affected. They run as separate processes with the Python named on their first line.

## Docker

The Docker images set this up already:

- `kaleido` comes from `docker/requirements.txt`, installed in the `python-base` image. An image built before kaleido was added does not have it. Rebuild and push `python-base` (`docker buildx bake --allow=fs.read=.. intermediate`), then rebuild `web`.
- Chromium is installed by the Playwright step in `docker/Dockerfile`.

See [Docker](./docker.md) for the build steps.

## Troubleshooting

The traceback for a failed download is in the Apache error log (`/var/log/apache2/error.log`, or the site's own error log, e.g. `ssl_umgear_error.log`). Match its last line:

| Error | Cause | Fix |
| --- | --- | --- |
| `ValueError: Image export using the "kaleido" engine requires the Kaleido package` | kaleido is not installed, or Apache was not restarted after installing it | [Step 1](#1-install-kaleido), then restart Apache |
| `ChromeNotFoundError` | No browser kaleido can find. Usually Chrome was installed with `plotly_get_chrome` without `--path`, into another user's home directory | [Step 2](#2-install-chrome-in-a-shared-location) |
| `BrowserFailedError: ... close immediately ...`, and the shell check in step 5 also fails | The browser cannot start: usually missing system libraries | Run `ldd` on the browser and install what is missing; then run it by hand as `www-data` with `--headless` (see [step 2](#2-install-chrome-in-a-shared-location)) to see its error |
| `BrowserFailedError: ... close immediately ...`, but the shell check in step 5 passes | The API's `sys.executable` is not the gEAR Python | Make sure the server has the current `www/api/api.wsgi`; see [Why api.wsgi sets sys.executable](#why-apiwsgi-sets-sysexecutable) |
| `BrowserFailedError` naming a browser you did not install in step 2 | kaleido picked a leftover browser | [Step 3](#3-check-which-browser-kaleido-uses) |

## Alternative: Playwright Chromium

[Playwright](https://playwright.dev/) publishes Chromium for both amd64 and arm64, and `playwright install --with-deps` also installs the system libraries it needs. The Docker image uses it (see `docker/Dockerfile`), because Docker on Apple-silicon Macs runs arm64, where Chrome for Testing is not available. Playwright itself is only needed for the download, so it goes in a throwaway virtual environment.

Use this instead of step 2, not as well as it: kaleido prefers Google Chrome, so any `/usr/bin/google-chrome` would be used instead (see [step 3](#3-check-which-browser-kaleido-uses)). Change `PYTHON` on the second line if your Python version differs:

```bash
sudo bash -euo pipefail -c '
PYTHON=/opt/Python-3.14.4/bin/python3
export PLAYWRIGHT_BROWSERS_PATH=/opt/ms-playwright
rm -rf /tmp/playwright-venv
$PYTHON -m venv /tmp/playwright-venv
/tmp/playwright-venv/bin/pip install -q playwright==1.55.0
/tmp/playwright-venv/bin/playwright install --with-deps --no-shell chromium
CHROME_BIN=$(find /opt/ms-playwright -type f -name chrome -path "*/chromium-*" | head -n 1)
[ -n "$CHROME_BIN" ] || { echo "ERROR: no Chromium under /opt/ms-playwright" >&2; exit 1; }
echo "Chromium: $CHROME_BIN"
ln -sf "$CHROME_BIN" /usr/bin/chromium
chmod -R a+rX /opt/ms-playwright
rm -rf /tmp/playwright-venv
/usr/bin/chromium --version'
```

The script stops at the first command that fails, and should end by printing a Chromium version. `--with-deps` supports Ubuntu and Debian releases only.

