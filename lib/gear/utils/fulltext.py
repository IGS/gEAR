"""
fulltext.py - Turn user search text into a MySQL full-text BOOLEAN MODE query.

In boolean mode, characters such as - ( ) + " are operators, so raw user text can exclude words
("Garcia-Añoveros" means "Garcia but not Añoveros") or be a syntax error ("Garcia-Añoveros)").
"""

import re

# A plain word (Unicode letters, digits, underscore), optionally with a trailing * wildcard
_PLAIN_WORD = re.compile(r"\w+\*?")
_NON_WORD = re.compile(r"[^\w]+")

# Shorter words aren't matched as prefixes or substrings: they would match most rows (e.g. "a*").
# Same as MySQL's default innodb_ft_min_token_size.
MIN_PARTIAL_MATCH_LENGTH = 3


def _with_prefix_wildcard(word: str) -> str:
    """Return word with a trailing * if it's long enough to match as a prefix."""
    if word.endswith("*") or len(word) < MIN_PARTIAL_MATCH_LENGTH:
        return word
    return word + "*"


def to_boolean_mode_query(text: str | None, prefix: bool = False) -> str:
    """
    Convert search text into a boolean-mode query that matches any of its terms.

    Whitespace-separated terms that are plain words (optionally ending in *) are kept as is.
    Any other term is matched as a phrase of its words, e.g. Garcia-Añoveros) -> "Garcia Añoveros".
    With prefix=True, single words of 3 or more characters also match as prefixes, so "Tono" finds
    "Tonotopic" (phrases can't take a wildcard). Returns "" if no searchable words remain.
    """
    terms = []
    for term in (text or "").split():
        if _PLAIN_WORD.fullmatch(term):
            terms.append(_with_prefix_wildcard(term) if prefix else term)
            continue
        words = _NON_WORD.sub(" ", term).split()
        if len(words) == 1:
            terms.append(_with_prefix_wildcard(words[0]) if prefix else words[0])
        elif words:
            terms.append('"{}"'.format(" ".join(words)))
    return " ".join(terms)


def to_like_patterns(text: str | None) -> list[str]:
    """
    Return a LIKE pattern ("%word%") for each word of 3 or more characters in the search text,
    without duplicates. Full-text search treats "_" as part of a word, so it can't find "potential" in
    "Groves_potential_targets"; these patterns can. LIKE's own wildcards (%, _ and the \\ escape) in
    the words are escaped.
    """
    patterns = []
    for word in _NON_WORD.sub(" ", text or "").split():
        if len(word) < MIN_PARTIAL_MATCH_LENGTH:
            continue
        escaped = word.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
        pattern = f"%{escaped}%"
        if pattern not in patterns:
            patterns.append(pattern)
    return patterns
