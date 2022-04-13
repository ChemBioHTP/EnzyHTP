"""Testing the constants and functions found in enzy_htp.chemical.residue.

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-19
"""
import pytest
from typing import List

from enzy_htp.core import InvalidResidueCode
from enzy_htp.chemical import residue as res

from util import all_caps


def test_variable_checks():
    """Basic checks on all of the globally defined variables in the enzy_htp.chemical.residue sub-module."""
    assert len(res.AA_LIST) == 21
    assert all_caps(res.AA_LIST)

    assert len(res.THREE_LETTER_AA_MAPPER) == 24

    three_mapper_keys = list(res.THREE_LETTER_AA_MAPPER.keys())
    three_mapper_values = list(set(list(res.THREE_LETTER_AA_MAPPER.values())))

    assert len(three_mapper_keys) == 24
    assert len(three_mapper_values) == 21

    assert all_caps(three_mapper_keys)
    assert all_caps(three_mapper_values)

    assert list(set(list(map(len, three_mapper_keys))))[0] == 3
    assert list(set(list(map(len, three_mapper_values))))[0] == 1

    assert len(res.ONE_LETTER_AA_MAPPER) == 21

    one_mapper_keys = list(res.ONE_LETTER_AA_MAPPER.keys())
    one_mapper_values = list(set(list(res.ONE_LETTER_AA_MAPPER.keys())))

    assert len(one_mapper_keys) == 21
    assert len(one_mapper_values) == 21

    assert all_caps(one_mapper_keys)
    assert all_caps(one_mapper_values)

    assert len(res.RESIDUE_ELEMENT_MAP) == 1
    assert "Amber" in res.RESIDUE_ELEMENT_MAP


def test_convert_to_three_letter():
    """Checking proper behavior and throwing of convert_to_three_letter()"""
    for one, three in res.ONE_LETTER_AA_MAPPER.items():
        assert res.convert_to_three_letter(one) == three
        assert res.convert_to_three_letter(one.lower()) == three

    with pytest.raises(InvalidResidueCode) as exe:
        res.convert_to_three_letter("")
    assert exe
    assert exe.type == InvalidResidueCode

    with pytest.raises(InvalidResidueCode) as exe:
        res.convert_to_three_letter("X")
    assert exe
    assert exe.type == InvalidResidueCode


def test_convert_to_one_letter():
    """Checking proper behavior and throwing of convert_to_one_letter()"""
    for three, one, in res.THREE_LETTER_AA_MAPPER.items():
        assert res.convert_to_one_letter(three) == one
        assert res.convert_to_one_letter(three.lower()) == one

    with pytest.raises(InvalidResidueCode) as exe:
        res.convert_to_one_letter("A")
    assert exe
    assert exe.type == InvalidResidueCode

    with pytest.raises(InvalidResidueCode) as exe:
        res.convert_to_one_letter("AXA")
    assert exe
    assert exe.type == InvalidResidueCode


def test_get_element_aliases():
    """Checking proper behavior and throwing of get_element_aliases()"""
    result = res.get_element_aliases("Amber", "H")
    assert result
    assert all_caps(result)

    with pytest.raises(SystemExit) as exe:
        res.get_element_aliases("DNE", "")

    assert exe
    assert exe.type == SystemExit

    assert not res.get_element_aliases("Amber", "ZDE")


def test_one_letters_except():
    """Checking that one_letters_except() corectly builds a list missing just the specified amino acid."""
    for aa in res.AA_LIST:
        local_res = res.one_letters_except(aa)
        assert len(local_res) == 20
        assert aa not in set(local_res)
        assert all_caps(local_res)

        local_res = res.one_letters_except(aa.lower())
        assert len(local_res) == 20
        assert aa not in set(local_res)
        assert all_caps(local_res)

    with pytest.raises(InvalidResidueCode) as exe:
        res.one_letters_except("X")

    assert exe
    assert exe.type == InvalidResidueCode
