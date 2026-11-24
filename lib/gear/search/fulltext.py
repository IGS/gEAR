"""
fulltext.py - Turn user search text into a MySQL full-text BOOLEAN MODE query.

In boolean mode, characters such as - ( ) + " are operators, so raw user text can exclude words
("Garcia-Añoveros" means "Garcia but not Añoveros") or be a syntax error ("Garcia-Añoveros)").

Searches match any of their terms. Text in double quotes is the one operator users can type: it
matches only as that exact phrase, so "hair cell" doesn't also match every dataset that mentions
"cell".
"""

import re

# A plain word (Unicode letters, digits, underscore), optionally with a trailing * wildcard
_PLAIN_WORD = re.compile(r"\w+\*?")
_NON_WORD = re.compile(r"[^\w]+")
# Text between a pair of double quotes. Curly quotes count too: macOS autocorrect and pasted text use them
_QUOTED = re.compile(r'["\u201c\u201d]([^"\u201c\u201d]*)["\u201c\u201d]')

# Shorter words aren't matched as prefixes or substrings: they would match most rows (e.g. "a*").
# Same as MySQL's default innodb_ft_min_token_size.
MIN_PARTIAL_MATCH_LENGTH = 3


def _with_prefix_wildcard(word: str) -> str:
    """Return word with a trailing * if it's long enough to match as a prefix."""
    if word.endswith("*") or len(word) < MIN_PARTIAL_MATCH_LENGTH:
        return word
    return word + "*"


def _split_on_quotes(text: str | None) -> list[tuple[bool, str]]:
    """
    Split search text into (is_quoted, chunk) pairs, in order. A lone, unbalanced quote isn't a
    pair, so it stays in an unquoted chunk and is later dropped as punctuation.
    """
    # re.split with one capture group alternates: unquoted, quoted, unquoted, ...
    chunks = _QUOTED.split(text or "")
    return [(index % 2 == 1, chunk) for index, chunk in enumerate(chunks)]


def _words(text: str) -> list[str]:
    """Return the words in text, with punctuation (including *) removed."""
    return _NON_WORD.sub(" ", text).split()


def _phrase(words: list[str]) -> str:
    return '"{}"'.format(" ".join(words))


def to_boolean_mode_query(text: str | None, prefix: bool = False, boost_phrase: bool = False) -> str:
    """
    Convert search text into a boolean-mode query that matches any of its terms.

    Whitespace-separated terms that are plain words (optionally ending in *) are kept as is.
    Any other term is matched as a phrase of its words, e.g. Garcia-Añoveros) -> "Garcia Añoveros".
    Text in double quotes matches exactly: several words as a phrase, one word as that whole word.

    With prefix=True, unquoted words of 3 or more characters also match as prefixes, so "Tono" finds
    "Tonotopic" (phrases can't take a wildcard).

    With boost_phrase=True, a search of two or more unquoted terms also ranks rows that contain all of
    them, in order, higher: hair cell -> hair* cell* >"hair cell". The same rows match either way;
    only their relevance score changes. Searches that use quotes are left as typed.

    Returns "" if no searchable words remain.
    """
    terms = []
    has_quotes = False
    unquoted_words = []    # every word of the unquoted terms, in order, for boost_phrase

    for is_quoted, chunk in _split_on_quotes(text):
        if is_quoted:
            words = _words(chunk)
            if words:
                has_quotes = True
                terms.append(_phrase(words) if len(words) > 1 else words[0])
            continue

        for term in chunk.split():
            if _PLAIN_WORD.fullmatch(term):
                terms.append(_with_prefix_wildcard(term) if prefix else term)
                unquoted_words.append(term.rstrip("*"))
                continue
            words = _words(term)
            if len(words) == 1:
                terms.append(_with_prefix_wildcard(words[0]) if prefix else words[0])
            elif words:
                terms.append(_phrase(words))
            unquoted_words.extend(words)

    if boost_phrase and not has_quotes and len(terms) > 1:
        terms.append(">" + _phrase(unquoted_words))
    return " ".join(terms)


def _escape_like(word: str) -> str:
    """Escape LIKE's own wildcards (%, _) and its \\ escape character."""
    return word.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def to_like_patterns(text: str | None) -> list[str]:
    """
    Return LIKE patterns that find the search text anywhere in a name, without duplicates.
    Full-text search treats "_" as part of a word, so it can't find "potential" in
    "Groves_potential_targets"; these patterns can. LIKE's own wildcards (%, _ and the \\ escape) in
    the words are escaped.

    Each unquoted word of 3 or more characters gets its own pattern ("%word%"). A quoted phrase gets
    one pattern with an unescaped _ (any one character) between its words, so "hair cell" also finds
    "Hair_Cell" and "hair-cell". A quoted single word gets no pattern, because quoting a word asks for
    that whole word, not for every name containing it.
    """
    patterns = []
    for is_quoted, chunk in _split_on_quotes(text):
        words = _words(chunk)
        if is_quoted:
            # A phrase of only very short words ("a b") would match most names
            if len(words) < 2 or max(len(word) for word in words) < MIN_PARTIAL_MATCH_LENGTH:
                continue
            candidates = ["%{}%".format("_".join(_escape_like(word) for word in words))]
        else:
            candidates = [f"%{_escape_like(word)}%" for word in words if len(word) >= MIN_PARTIAL_MATCH_LENGTH]

        for pattern in candidates:
            if pattern not in patterns:
                patterns.append(pattern)
    return patterns
