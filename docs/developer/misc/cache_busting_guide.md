# Cache Busting Implementation Guide

## Overview

gEAR prevents stale browser caches by appending a version query string (`?v=<cache_version>`) to page-specific CSS and JavaScript URLs. The version is bumped automatically on every git commit.

## How It Works

1. **Version storage**: The version lives in `www/cache_version.json`:

   ```json
   {
       "cache_version": "2026.09.03.212437"
   }
   ```

2. **Auto-bump**: A pre-commit hook (`.githooks/pre-commit`) replaces the value with a timestamp whenever other files are staged.
3. **Loading**: `getDomainPreferences()` in `www/js/common.v2.js` fetches `/site_domain_prefs.json` and `/cache_version.json` and merges the version into the returned object as `prefs.cache_version`. Pages then load their assets with `insertVersionedCSS()` / `insertVersionedJS()`.

## Helpers in `www/js/common.v2.js`

| Function | Behavior |
| --- | --- |
| `getDomainPreferences()` | Returns site preferences with `cache_version` merged in |
| `insertVersionedCSS(href, cacheVersion)` | Appends `<link rel="stylesheet" href="<href>?v=<cacheVersion>">` to `<head>` |
| `insertVersionedJS(href, cacheVersion)` | Appends `<script type="module" src="<href>?v=<cacheVersion>">` to `<body>` |

`insertVersionedJS` always creates a module script.

## Usage

Nearly every page in `www/` ends with an inline module like this one from `www/dataset_curator.html`:

```html
<!-- Page-specific CSS/JS loading here -->
<script type="module">
  import { getDomainPreferences, insertVersionedCSS, insertVersionedJS, loadDomainFunding } from './js/common.v2.js';
  const prefs = await getDomainPreferences();

  await loadDomainFunding(prefs);

  // Load any local JS with cache-busting version parameter
  insertVersionedJS('js/dataset_curator.js', prefs.cache_version);

  // Now load page-specific CSS
  insertVersionedCSS('css/common.v2.css', prefs.cache_version);
  insertVersionedCSS('css/curator_common.css', prefs.cache_version);
  insertVersionedCSS('css/dataset_curator.css', prefs.cache_version);
</script>
```

When adding a new page, copy this block and replace the page-specific file names. Vendor libraries and CDN assets are loaded with ordinary `<script>`/`<link>` tags and are not versioned. `common.v2.js` itself is imported without a version string.

Pages that currently do not use the helpers: the legacy `contact.html`.

## Git Hook

The hook is registered with the [pre-commit](https://pre-commit.com/) framework in `.pre-commit-config.yaml` (hook id `cache-version-bump`, entry `.githooks/pre-commit`). Install it once per clone:

```bash
pip install pre-commit
pre-commit install
```

On each commit the hook:

- Generates a timestamp version (`%Y.%m.%d.%H%M%S`, e.g. `2026.02.10.141530`)
- Updates `www/cache_version.json` and stages it
- Skips the bump if `cache_version.json` is the only staged file (prevents loops)
- Skips the bump if `cache_version.json` has unstaged edits

See [Developer Documentation](../README.md#setting-up-a-development-environment) for general setup.

### Testing the Hook

```bash
echo "test" >> test.txt
git add test.txt
git commit -m "Test cache version bump"
# ✓ Bumped cache_version: 2026.02.10.141530 → 2026.02.10.141812
```

## Auditing

`bin/audit_cache_busting.py` scans `www/**/*.html` for local CSS/JS references without `?v=`. It reads the current version from `www/cache_version.json` and counts a page as compliant if it calls `insertVersionedJS()` or `insertVersionedCSS()`. Partials under `www/include/` are loaded into pages rather than served directly, so they show up in the "needs attention" list even though they don't need their own versioning.
