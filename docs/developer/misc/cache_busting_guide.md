# Cache Busting Implementation Guide

## Overview
This system automatically prevents browser caching issues by appending a version number to CSS and JavaScript asset URLs. The version is automatically updated on every git commit.

## How It Works

1. **Version Storage**: The cache version is stored in `www/site_domain_prefs.json` as `cache_version`
2. **Auto-Bump**: A git pre-commit hook automatically updates the version timestamp whenever you commit (only when there are code changes, not for cache_version-only commits)
3. **Dynamic Loading**: JavaScript can append this version to asset URLs to force cache refresh

## Usage Examples

### Option 1: Dynamic JS Loading (Recommended for JS modules)

For pages that use ES modules, you can dynamically load and version assets:

```javascript
import { getDomainPreferences, versionedAsset } from './js/common.v2.js';

// Load site preferences first
const prefs = await getDomainPreferences();

// Dynamically load a CSS file with versioning
const link = document.createElement('link');
link.rel = 'stylesheet';
link.href = versionedAsset('css/my-styles.css');
document.head.appendChild(link);

// Dynamically load a JS file with versioning
const script = document.createElement('script');
script.src = versionedAsset('js/my-script.js');
document.body.appendChild(script);
```

### Option 2: Template Variable (Recommended for HTML)

If you're using a template engine (Jinja2, Mako, etc.) in your Python backend:

**In your Python view/controller:**
```python
import json

def load_cache_version():
    with open('www/site_domain_prefs.json') as f:
        prefs = json.load(f)
        return prefs.get('cache_version', '1.0.0')

# Pass to template
return render_template('index.html', cache_version=load_cache_version())
```

**In your HTML template:**
```html
<link rel="stylesheet" href="css/common.v2.css?v={{ cache_version }}" />
<script type="module" src="js/index.js?v={{ cache_version }}"></script>
```

### Option 3: Inline JavaScript (For static HTML)


## Git Hook

The pre-commit hook (defined in `.githooks/pre-commit`) automatically:
- Generates a new timestamp-based version (e.g., `2026.02.10.141530`)
- Updates `www/site_domain_prefs.json`
- Stages the updated file in your commit
- **Smart mode**: Skips the update if `site_domain_prefs.json` is the only staged change, preventing infinite loops

The hook is managed through the [pre-commit](https://pre-commit.com/) framework. See [Developer Documentation](./developer/README.md#setting-up-a-development-environment) for setup instructions.

### Testing the Hook

```bash
# Make a change and commit
echo "test" >> test.txt
git add test.txt
git commit -m "Test cache version bump"
# You should see: ✓ Bumped cache_version: 1.0.0 → 2026.02.10.141530

# On the follow-up commit of just the cache_version change, you'll see:
# ℹ Skipping cache_version bump (only cache_version changed)
```
