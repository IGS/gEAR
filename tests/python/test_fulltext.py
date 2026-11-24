"""gear.utils.fulltext: converting search text into a MySQL boolean-mode query."""

import pytest

from gear.utils.fulltext import to_boolean_mode_query, to_like_patterns


@pytest.mark.parametrize("text, expected", [
    # "-" would exclude Añoveros, and an unbalanced ")" is a syntax error
    ("Garcia-Añoveros", '"Garcia Añoveros"'),
    ("Garcia-Añoveros)", '"Garcia Añoveros"'),
    ("Garcia", "Garcia"),
    ("Añoveros", "Añoveros"),
    ("hair cell", "hair cell"),
    ("  hair   cell ", "hair cell"),
    ("cochlea*", "cochlea*"),
    ("E18.5-P0", '"E18 5 P0"'),
    ("GSE137299", "GSE137299"),
    ("-foo +bar", "foo bar"),
    ('"hair cell"', '"hair cell"'),                  # quoted text is an exact phrase
    ("(INSM1) @3 ~x <y >z", "INSM1 3 x y z"),
    ("co*chlea", '"co chlea"'),
    ("Müller glia", "Müller glia"),
])
def test_conversion(text, expected):
    assert to_boolean_mode_query(text) == expected


@pytest.mark.parametrize("text", [None, "", "   ", "()", "-", '"', '""', "* + -"])
def test_nothing_searchable(text):
    assert to_boolean_mode_query(text) == ""


@pytest.mark.parametrize("text, expected", [
    # #636: words of 3+ characters match as prefixes
    ("Tono", "Tono*"),
    ("hair cell", "hair* cell*"),
    ("SURFACE", "SURFACE*"),
    ("Groves_potential", "Groves_potential*"),
    ("cochlea*", "cochlea*"),
    ("ab P0", "ab P0"),                              # too short to match as a prefix
    ("Garcia-Añoveros", '"Garcia Añoveros"'),        # a phrase can't take a wildcard
    ("(INSM1)", "INSM1*"),                            # one word left after removing punctuation
])
def test_prefix_conversion(text, expected):
    assert to_boolean_mode_query(text, prefix=True) == expected


def test_prefix_is_off_by_default():
    assert to_boolean_mode_query("Tono") == "Tono"


@pytest.mark.parametrize("text, expected", [
    # Quoted text matches exactly: as a phrase, or as a whole word with no prefix wildcard
    ('"hair cell"', '"hair cell"'),
    ('"cell"', "cell"),
    ('"hair cell" cochlea', '"hair cell" cochlea*'),
    ('cochlea "hair cell"', 'cochlea* "hair cell"'),
    ("\u201chair cell\u201d", '"hair cell"'),             # curly quotes, e.g. from macOS autocorrect
    ('"Garcia-Añoveros"', '"Garcia Añoveros"'),
    ('"hair cell', "hair* cell*"),                      # an unbalanced quote is ignored
    ('"hair cell" ""', '"hair cell"'),                  # empty quotes add nothing
])
def test_quoted_phrases(text, expected):
    assert to_boolean_mode_query(text, prefix=True) == expected


@pytest.mark.parametrize("text, expected", [
    # Several unquoted terms also rank rows containing all of them, in order, higher (">")
    ("hair cell", 'hair* cell* >"hair cell"'),
    ("outer hair cell", 'outer* hair* cell* >"outer hair cell"'),
    ("Garcia-Añoveros cochlea", '"Garcia Añoveros" cochlea* >"Garcia Añoveros cochlea"'),
    ("Tono", "Tono*"),                                  # a single word has no phrase to boost
    ("Garcia-Añoveros", '"Garcia Añoveros"'),           # already a phrase
    ('"hair cell" cochlea', '"hair cell" cochlea*'),    # the user chose a phrase with quotes
])
def test_phrase_boost(text, expected):
    assert to_boolean_mode_query(text, prefix=True, boost_phrase=True) == expected


def test_phrase_boost_is_off_by_default():
    assert to_boolean_mode_query("hair cell", prefix=True) == "hair* cell*"


@pytest.mark.parametrize("text, expected", [
    ("potential", ["%potential%"]),
    ("Groves_potential targets", ["%Groves\\_potential%", "%targets%"]),   # _ is a LIKE wildcard
    ("Garcia-Añoveros", ["%Garcia%", "%Añoveros%"]),
    ("ab P0 x", []),                                                        # too short
    ("hair hair*", ["%hair%"]),                                             # no duplicates
    ("50% 100\\x", ["%100%"]),                                             # % and \ never reach a pattern
    # A quoted phrase is one pattern; the unescaped _ between words also matches "Hair_Cell" or "hair-cell"
    ('"hair cell"', ["%hair_cell%"]),
    ('"hair cell" cochlea', ["%hair_cell%", "%cochlea%"]),
    ('"Groves_potential targets"', ["%Groves\\_potential_targets%"]),
    ('"cell"', []),                                                         # one quoted word: whole word only
    ('"a b"', []),                                                          # would match most names
    (None, []),
    ("", []),
])
def test_like_patterns(text, expected):
    assert to_like_patterns(text) == expected
