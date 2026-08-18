#!/opt/bin/python3
'''
Reverse of the www/p script: takes a long-form gEAR URL and returns the /p shorthand URL.

Usage:
    python3 create_shortform_url.py --url "https://umgear.org/expression.html?share_id=abc123&gene_symbol=Atoh1"

To validate the roundtrip, pass the returned /p URL back into the /p script
and compare the resulting destination URL with your original input.

Long-form parameter mappings:
    share_id                -> s
    layout_id               -> l
    gene_symbol             -> g
    gene_symbol_exact_match -> gsem
    gene_lists              -> c  (non-projection pages)
    projection_source       -> c  (projection page)
    is_multigene            -> multi  (non-projection pages)
    multipattern_plots      -> multi  (projection page)
    projection_patterns     -> ptrns
    projection_algorithm    -> algo
    zscore                  -> zscore
    minclip                 -> minclip
    search_string           -> ss  (de/gcm pages only)

Parameters not mapped:
    identifiers_org_id (id) - This is meant for unpacking identifiers.org URIs.
    page (p) - This is implied by the page filename in the long-form URL.

'''

import argparse
import urllib.parse


# Map long-form page filenames/names to /p 'p' parameter values
PAGE_NAME_TO_SHORT = {
    'expression.html': None,        # default; omit 'p' param
    'sc_workbench.html': 'a',
    'dataset_explorer.html': 'de',
    'gene_list_manager.html': 'gcm',
    'projection.html': 'p',
}


def check_for_truthy(value: str | None) -> int:
    """Return 1 if value is truthy, 0 otherwise (mirrors logic in www/p)."""
    if value is None:
        return 0
    if isinstance(value, str) and value.lower() in ('0', 'false', 'f', 'no', 'n'):
        return 0
    return 1


def resolve_page_short(page_name: str) -> str | None:
    """
    Resolve the short page code from a full or bare page name.

    Args:
        page_name: e.g. 'expression.html', 'expression', 'sc_workbench', etc.

    Returns:
        Short code string (e.g. 'a', 'de') or None for the default expression page.
    """
    normalised = page_name.lower().strip()
    for full_name, short in PAGE_NAME_TO_SHORT.items():
        if normalised in (full_name, full_name.replace('.html', '')):
            return short
    # Unrecognised page: default to expression
    return None


def build_short_params(query_params: dict, page_short: str | None) -> dict:
    """
    Convert long-form query parameters into /p short-form parameters.

    Args:
        query_params: Dict of long-form param name -> value (first value only).
        page_short: Short page identifier (None for expression.html default).

    Returns:
        Ordered dict of short-form parameter key/value pairs.
    """
    params = {}

    if page_short is not None:
        params['p'] = page_short

    share_id = query_params.get('share_id')
    if share_id is not None:
        params['s'] = share_id

    layout_id = query_params.get('layout_id')
    if layout_id is not None:
        params['l'] = layout_id

    gene_symbol = query_params.get('gene_symbol')
    if gene_symbol is not None:
        params['g'] = gene_symbol

    # gsem only applies to non-special pages (mirrors www/p logic)
    if page_short not in ('a', 'de', 'gcm', 'p'):
        gsem = query_params.get('gene_symbol_exact_match', '1')
        params['gsem'] = check_for_truthy(gsem)

    # gene cart: projection uses 'projection_source', others use 'gene_lists'
    gene_cart_key = 'projection_source' if page_short == 'p' else 'gene_lists'
    gene_cart = query_params.get(gene_cart_key)
    if gene_cart is not None:
        params['c'] = gene_cart

    # multigene: projection uses 'multipattern_plots', others use 'is_multigene'
    multi_key = 'multipattern_plots' if page_short == 'p' else 'is_multigene'
    multi = query_params.get(multi_key)
    if multi is not None:
        params['multi'] = multi

    # Projection-specific params
    if page_short == 'p':
        patterns = query_params.get('projection_patterns')
        if patterns is not None:
            params['ptrns'] = patterns

        algo = query_params.get('projection_algorithm')
        if algo is not None:
            params['algo'] = algo

        params['zscore'] = check_for_truthy(query_params.get('zscore', '0'))

        minclip = query_params.get('minclip')
        if minclip is not None:
            params['minclip'] = minclip

    # search_string only for de / gcm pages
    if page_short in ('de', 'gcm'):
        search_string = query_params.get('search_string')
        if search_string is not None:
            params['ss'] = urllib.parse.quote(search_string)

    return params


def convert_url(long_url: str) -> str:
    """
    Convert a long-form gEAR URL to its /p shorthand equivalent.

    Args:
        long_url: Full gEAR URL, e.g.
            'https://umgear.org/expression.html?share_id=abc&gene_symbol=Atoh1'

    Returns:
        The equivalent /p shorthand URL string.
    """
    parsed = urllib.parse.urlparse(long_url)

    # Derive the base URL up to (but not including) the page filename
    # e.g. https://umgear.org/expression.html -> https://umgear.org/
    path = parsed.path.rstrip('/')
    base = '{scheme}://{host}{pre}/p'.format(
        scheme=parsed.scheme or 'https',
        host=parsed.netloc,
        pre=path.rsplit('/', 1)[0] if '/' in path else '',
    )

    # Identify the page from the path filename
    page_filename = path.rsplit('/', 1)[-1]  # e.g. 'expression.html'
    page_short = resolve_page_short(page_filename)

    # Parse query string into a flat dict (first value wins for duplicates)
    query_params = {
        k: v[0]
        for k, v in urllib.parse.parse_qs(parsed.query, keep_blank_values=True).items()
    }

    short_params = build_short_params(query_params, page_short)

    if short_params:
        query_string = urllib.parse.urlencode(short_params)
        return '{0}?{1}'.format(base, query_string)

    return base


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Convert a long-form gEAR URL to its /p shorthand equivalent.',
    )
    parser.add_argument(
        '--url',
        required=True,
        metavar='URL',
        help=(
            'Full long-form gEAR URL to convert, e.g. '
            '"https://umgear.org/expression.html?share_id=abc&gene_symbol=Atoh1"'
        ),
    )
    args = parser.parse_args()

    short_url = convert_url(args.url)
    print(short_url)


if __name__ == '__main__':
    main()