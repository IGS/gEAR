"""
fulltext.py - Turn user search text into a MySQL full-text BOOLEAN MODE query.

In boolean mode, characters such as - ( ) + " are operators, so raw user text can exclude words
("Garcia-Añoveros" means "Garcia but not Añoveros") or be a syntax error ("Garcia-Añoveros)").
"""

import re

# A plain word (Unicode letters, digits, underscore), optionally with a trailing * wildcard
_PLAIN_WORD = re.compile(r"\w+\*?")
_NON_WORD = re.compile(r"[^\w]+")


def to_boolean_mode_query(text: str | None) -> str:
    """
    Convert search text into a boolean-mode query that matches any of its terms.

    Whitespace-separated terms that are plain words (optionally ending in *) are kept as is.
    Any other term is matched as a phrase of its words, e.g. Garcia-Añoveros) -> "Garcia Añoveros".
    Returns "" if no searchable words remain.
    """
    terms = []
    for term in (text or "").split():
        if _PLAIN_WORD.fullmatch(term):
            terms.append(term)
            continue
        words = _NON_WORD.sub(" ", term).split()
        if len(words) == 1:
            terms.append(words[0])
        elif words:
            terms.append('"{}"'.format(" ".join(words)))
    return " ".join(terms)
