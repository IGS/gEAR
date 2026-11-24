"""
Keyword matching in Python, for searchers whose text isn't in MySQL (study pages, cell-cell
communication analyses).

The rules match the MySQL full-text search used for datasets and gene lists
(gear.utils.fulltext), so a search behaves the same for every type:
- any one term is enough to match,
- unquoted words of 3 or more characters also match as prefixes ("cochle" finds "cochlear"),
- quoted text matches only as that phrase, or as that whole word.
"""

import re

# The quote splitting and word rules are shared with the SQL search, so both read
#  search text the same way
from gear.utils.fulltext import MIN_PARTIAL_MATCH_LENGTH, _split_on_quotes, _words


def compile_terms(text):
    """Return one case-insensitive regex per search term. An empty list means no keywords."""
    patterns = []
    for is_quoted, chunk in _split_on_quotes(text):
        words = _words(chunk)
        if is_quoted:
            if words:
                phrase = r"\W+".join(re.escape(word) for word in words)
                patterns.append(re.compile(rf"\b{phrase}\b", re.IGNORECASE))
            continue

        for word in words:
            # Short words match only as whole words, or they would match almost everything
            end = "" if len(word) >= MIN_PARTIAL_MATCH_LENGTH else r"\b"
            patterns.append(re.compile(rf"\b{re.escape(word)}{end}", re.IGNORECASE))
    return patterns


def score(patterns, weighted_fields):
    """
    Score text against the search terms. 0 means no term matched.

    weighted_fields is a list of (text, weight). Each term adds the weight of the most important
    field it appears in, so a title match counts for more than a match in the body text.
    """
    total = 0
    for pattern in patterns:
        total += max((weight for text, weight in weighted_fields if text and pattern.search(text)), default=0)
    return total
