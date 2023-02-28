"""Testing the enzy_htp.mutation.mutation_pattern.target_aa_pattern.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-17
"""

import os
import enzy_htp.mutation.mutation_pattern.target_aa_pattern as m_p

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"


def test_decode_target_aa_pattern():
    """test the function use a made up target_aa_pattern"""
    test_pattern = "self"
    test_orig_resi = "ARG"
    assert m_p.decode_target_aa_pattern(test_orig_resi, test_pattern) == set(["ARG"])


def test_decode_target_aa_pattern_not_expr():
    """test the function use a made up target_aa_pattern"""
    test_pattern = "all not self"
    test_orig_resi = "ARG"
    assert m_p.decode_target_aa_pattern(test_orig_resi, test_pattern) == set([
        "GLU", "ASN", "HIS", "SER", "VAL", "GLN", "TYR", "PRO", "ASP", "CYS", "GLY",
        "ALA", "THR", "LYS", "LEU", "MET", "ILE", "TRP", "PHE"
    ])


def test_decode_larger():
    """test the function use a made up orig_resi"""
    test_orig_resi = "ARG"
    assert m_p.decode_larger(test_orig_resi) == ['PHE', 'TRP', 'TYR']


def test_decode_charge_p1():
    """test the function use a made up orig_resi"""
    test_orig_resi = "GLY"
    assert m_p.decode_charge_n1(test_orig_resi) == ['ASP', 'GLU']


def test_check_target_aa_pattern_charge(caplog):
    """test the function give correct warning about HIS"""
    test_pattern = "charge+"
    m_p.check_target_aa_pattern(test_pattern)
    assert "HIS is not included" in caplog.text
