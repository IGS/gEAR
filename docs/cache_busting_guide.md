# Cache Busting Implementation Guide

## Overview
This system automatically prevents browser caching issues by appending a version number to CSS and JavaScript asset URLs. The version is automatically updated on every git commit.

## How It Works

1. **Version Storage**: The cache version is stored in `www/site_domain_prefs.json` as `cache_version`
2. **Auto-Bump**: A git pre-commit hook automatically updates the version timestamp whenever you commit
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

Add this script early in your HTML `<head>`:

```html
<script>
  // Load cache version and update asset URLs
  (async function() {
    const response = await fetch('/site_domain_prefs.json');
    const prefs = await response.json();
    const version = prefs.cache_version;
    
    // Update all local stylesheet hrefs
    document.querySelectorAll('link[rel="stylesheet"]').forEach(link => {
      const href = link.getAttribute('href');
      if (href && !href.startsWith('http') && !href.includes('?v=')) {
        link.href = `${href}?v=${version}`;
      }
    });
    
    // Update all local script sources
    document.querySelectorAll('script[src]').forEach(script => {
      const src = script.getAttribute('src');
      if (src && !src.startsWith('http') && !src.includes('?v=')) {
        script.src = `${src}?v=${version}`;
      }
    });
  })();
</script>

<!-- Your regular asset links -->
<link rel="stylesheet" href="css/common.v2.css" />
<script src="js/vendor/js.cookie.js"></script>
```

### Option 4: Build-time Replacement

Create a deployment script that replaces a placeholder:

**In your HTML:**
```html
<link rel="stylesheet" href="css/common.v2.css?v=__CACHE_VERSION__" />
<script src="js/index.js?v=__CACHE_VERSION__"></script>
```

**Deployment script:**
```bash
#!/bin/bash
VERSION=$(grep -Po '"cache_version":\s*"\K[^"]*' www/site_domain_prefs.json)

# Replace placeholders in HTML files
find www/ -name "*.html" -exec sed -i "s/__CACHE_VERSION__/${VERSION}/g" {} \;
```

## Git Hook

The pre-commit hook (`.git/hooks/pre-commit`) automatically:
- Generates a new timestamp-based version (e.g., `2026.02.10.141530`)
- Updates `www/site_domain_prefs.json`
- Stages the updated file in your commit

### Testing the Hook

```bash
# Make a change and commit
echo "test" >> test.txt
git add test.txt
git commit -m "Test cache version bump"
# You should see: ✓ Bumped cache_version: 1.0.0 → 2026.02.10.141530
```

## Best Practices

1. **Don't version external CDN assets** - They manage their own versioning
2. **Version all local CSS/JS files** - Ensures users get updates
3. **HTML files should NOT be cached** - Use proper cache headers in Apache
4. **Test after deployment** - Verify version parameter appears in browser DevTools

## Apache Configuration

Add to your `.htaccess` or Apache config:

```apache
# Never cache HTML
<FilesMatch "\.(html|htm)$">
    Header set Cache-Control "no-cache, no-store, must-revalidate"
    Header set Pragma "no-cache"
    Header set Expires "0"
</FilesMatch>

# Cache versioned assets for 1 year
<FilesMatch "\.(css|js|jpg|png|gif|svg|woff|woff2)$">
    Header set Cache-Control "max-age=31536000, public"
</FilesMatch>
```

## Troubleshooting

**Version not updating?**
- Check git hook is executable: `ls -l .git/hooks/pre-commit`
- Check hook output: Look for "✓ Bumped cache_version" message

**Assets still cached?**
- Hard refresh browser: Ctrl+F5 (Windows/Linux) or Cmd+Shift+R (Mac)
- Check DevTools Network tab - version parameter should be in URL
- Verify `site_domain_prefs.json` has current version

**Hook not running?**
- Ensure you're committing changes: `git add . && git commit`
- Check hook has execute permissions: `chmod +x .git/hooks/pre-commit`
# Test commit to trigger version bump
