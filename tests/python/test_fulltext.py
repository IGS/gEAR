"""gear.utils.fulltext: converting search text into a MySQL boolean-mode query."""

import pytest

from gear.utils.fulltext import to_boolean_mode_query


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
    ('"hair cell"', 'hair cell'),
    ("(INSM1) @3 ~x <y >z", "INSM1 3 x y z"),
    ("co*chlea", '"co chlea"'),
    ("Müller glia", "Müller glia"),
])
def test_conversion(text, expected):
    assert to_boolean_mode_query(text) == expected


@pytest.mark.parametrize("text", [None, "", "   ", "()", "-", '"', "* + -"])
def test_nothing_searchable(text):
    assert to_boolean_mode_query(text) == ""
