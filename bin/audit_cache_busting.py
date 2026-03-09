#!/usr/bin/env python3
"""
Audit HTML files for cache-busting compliance.
Checks that local CSS/JS assets have versioning applied.
"""

import re
import json
from pathlib import Path

def load_cache_version():
    """Load current cache version from site_domain_prefs.json"""
    prefs_file = Path(__file__).parent.parent / 'www' / 'site_domain_prefs.json'
    with open(prefs_file) as f:
        prefs = json.load(f)
        return prefs.get('cache_version', 'UNKNOWN')

def check_html_file(html_path):
    """Check if HTML file properly versions its assets"""
    with open(html_path) as f:
        content = f.read()
    
    issues = []
    
    # Check for local CSS files without version
    css_pattern = r'<link[^>]*href=["\'](?!https?://)([^"\']+\.css)(?!\?v=)["\']'
    unversioned_css = re.findall(css_pattern, content)
    
    # Check for local JS files without version (excluding type="module" which may use dynamic loading)
    js_pattern = r'<script[^>]*src=["\'](?!https?://)([^"\']+\.js)(?!\?v=)["\'][^>]*(?!type=["\']module)'
    unversioned_js = re.findall(js_pattern, content)
    
    # Check if file has the inline versioning script
    has_inline_script = 'fetch(\'/site_domain_prefs.json\')' in content
    
    # Check if file uses versionedAsset function
    uses_versioned_asset = 'versionedAsset(' in content
    
    return {
        'unversioned_css': unversioned_css,
        'unversioned_js': unversioned_js,
        'has_inline_script': has_inline_script,
        'uses_versioned_asset': uses_versioned_asset,
        'has_cache_busting': has_inline_script or uses_versioned_asset
    }

def main():
    """Main audit function"""
    www_dir = Path(__file__).parent.parent / 'www'
    html_files = list(www_dir.rglob('*.html'))
    
    print(f"🔍 Auditing {len(html_files)} HTML files for cache-busting compliance")
    print(f"📦 Current cache version: {load_cache_version()}\n")
    
    compliant = []
    non_compliant = []
    
    for html_file in html_files:
        rel_path = html_file.relative_to(www_dir)
        result = check_html_file(html_file)
        
        if result['has_cache_busting']:
            compliant.append(rel_path)
        elif result['unversioned_css'] or result['unversioned_js']:
            non_compliant.append((rel_path, result))
    
    # Print results
    print("✅ Compliant files (using cache-busting):")
    for path in sorted(compliant):
        print(f"  ✓ {path}")
    
    print(f"\n⚠️  Files needing attention ({len(non_compliant)}):")
    for path, result in sorted(non_compliant):
        print(f"  ⚠️  {path}")
        if result['unversioned_css']:
            print(f"      Unversioned CSS: {', '.join(result['unversioned_css'][:3])}")
        if result['unversioned_js']:
            print(f"      Unversioned JS: {', '.join(result['unversioned_js'][:3])}")
        print(f"      Suggested fix: Add inline versioning script or use versionedAsset()")
    
    print(f"\n📊 Summary:")
    print(f"  Total files: {len(html_files)}")
    print(f"  Compliant: {len(compliant)} ({len(compliant)/len(html_files)*100:.1f}%)")
    print(f"  Need attention: {len(non_compliant)} ({len(non_compliant)/len(html_files)*100:.1f}%)")
    
    if non_compliant:
        print(f"\n💡 Tip: See docs/developer/cache_busting_guide.md for implementation options")
        return 1
    else:
        print(f"\n🎉 All HTML files are using cache-busting!")
        return 0

if __name__ == '__main__':
    exit(main())
