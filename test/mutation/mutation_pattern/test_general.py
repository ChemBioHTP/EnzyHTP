"""Testing the enzy_htp.mutation.mutation_pattern.general.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-01-26
"""

import numpy as np
import pytest
import os
import pickle

from enzy_htp.core.exception import InvalidMutationPatternSyntax
from enzy_htp import PDBParser
from enzy_htp.mutation.mutation import Mutation
import enzy_htp.mutation.mutation_pattern.general as m_p

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

@pytest.mark.interface
def test_decode_mutation_pattern():
    """dev run of the function"""
    test_mutation_pattern = (
        "KA162A, {RA154W, HA201A},"
        " r:2[resi 289 around 4 and not resi 36:larger,"
            " proj(ID 1000, ID 2023, positive, 10):more_negative_charge]*100"
        )
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)

    m_p.decode_mutation_pattern(test_stru, test_mutation_pattern)

def test_seperate_mutants():
    """test the function use a made up mutation_pattern for KE
    use the len of the seperation result as the fingerprint to assert"""
    test_mutation_pattern = (
        "KA162A, {RA154W, HA201A},"
        " r:2[resi 289 around 4 and not resi 36:larger,"
            " proj(ID 1000, ID 2023, positive, 10):more_negative_charge]*100"
        )
    assert len(m_p.seperate_mutant_patterns(test_mutation_pattern)) == 3

def test_seperate_sections():
    """test the function use a made up mutation_pattern for KE
    use the len of the seperation result as the fingerprint to assert"""
    test_mutation_pattern = (
        "KA162A,"
        " r:2[resi 289 around 4 and not resi 36:larger,"
            " proj(ID 1000, ID 2023, positive, 10):more_negative_charge]*100,"
        "RA154W"
        )
    assert len(m_p.seperate_section_patterns(test_mutation_pattern)) == 3

def test_get_section_type():
    """test if the section type is correctly determined"""
    section_1 = "KA162A"
    section_2 = "r:2[resi 289 around 4 and not resi 36:larger]"
    section_3 = "a:[resi 289 around 4 and not resi 36:larger]"

    assert m_p.get_section_type(section_1) == "d"
    assert m_p.get_section_type(section_2) == "r"
    assert m_p.get_section_type(section_3) == "a"

def test_get_section_type_bad():
    """test if the section type raise the execption as expected."""
    section_bad = "x:KA162A"
    with pytest.raises(Exception) as exe:
        m_p.get_section_type(section_bad)
    assert exe.type == InvalidMutationPatternSyntax

def test_decode_direct_mutation():
    """test the function works as expected using a made up pattern and manually
    curated answer. test giving default chain id"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_d_pattern = "RA154W"
    test_d_pattern_1 = "R154W"
    answer = Mutation(orig="R", target="W", chain_id="A", res_idx=154)

    assert m_p.decode_direct_mutation(test_stru, test_d_pattern) == answer
    assert m_p.decode_direct_mutation(test_stru, test_d_pattern_1) == answer

def test_decode_mutation_esm_pattern():
    """test function works as expected, test using pickle obj of a confirmed run"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 254 around 5:all not self, resi 2:larger"
    answer_path = f"{DATA_DIR}mutation_esm_answer.pickle"
    mutation_esm = m_p.decode_mutation_esm_pattern(test_stru, test_pattern)
    # with open(answer_path, "wb") as f:
    #     pickle.dump(mutation_esm, f, 0)
    with open(answer_path, "rb") as f:
        answer = pickle.load(f)
    assert answer == mutation_esm

def test_decode_mutation_esm_pattern_share_point():
    """test the case that there are shared positions from different
    mutation_esm_pattern s"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 9:all not self, resi 254 around 3:charge+"
    mutation_esm = m_p.decode_mutation_esm_pattern(test_stru, test_pattern)
    assert len(mutation_esm[('A', 9)]) == 19

def test_decode_random_mutation(caplog):
    """test the function works as expected using a made up pattern and manually
    curated answer. test of non repeat case
    Use a random seed to control the test to contain a repeating random result"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "r:2[resi 254 around 3:all not self]*10"
    # find a seed that have repeating mutant
    # i = True
    # while i:
    #     seed = np.random.randint(1, 10000)
    #     print(f"seed: {seed}")
    #     np.random.seed(seed)
    #     mutants = m_p.decode_random_mutation(test_stru, test_pattern)
    #     if "repeating" in caplog.text:
    #         i=False
    np.random.seed(457)
    mutants = m_p.decode_random_mutation(test_stru, test_pattern)
    assert len(mutants) == 10
    for i in mutants:
        assert len(i) == 2
    assert "repeating MUTANT is generated" in caplog.text

def test_decode_random_mutation_allow_repeat(caplog):
    """test the function works as expected using a made up pattern and manually
    curated answer. test of non repeat case
    Use a random seed to control the test to contain a repeating random result"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "r:2R[resi 254 around 3:all not self]*10R"
    np.random.seed(457) # seed that contains a repeating mutant

    mutants = m_p.decode_random_mutation(test_stru, test_pattern)
    assert len(mutants) == 10
    for i in mutants:
        assert len(i) > 0
    assert "repeating mutation is generated" in caplog.text

def test_decode_all_mutation():
    """test the function works as expected using a made up pattern and manually
    curated answer."""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "a:[resi 253:all not self, resi 252:larger]"

    mutants = m_p.decode_all_mutation(test_stru, test_pattern)
    assert len(mutants) == 400

def test_decode_all_mutation_m_flag():
    """test the function with the flag M specificed
    works as expected using a made up pattern and manually
    curated answer."""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "a:M[resi 253:all not self, resi 252:larger]"

    mutants = m_p.decode_all_mutation(test_stru, test_pattern)
    assert len(mutants) == 361

def test_combine_section_mutant_one_to_many():
    """test the function works as expected in the case that
    single mutant combine with many mutants"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_sec_1 = "a:[resi 253:all not self, resi 252:larger]"
    test_sec_3 = "L10A"
    # build per-section mutant mapper
    p_mutant_mapper = {}
    p_mutant_mapper[test_sec_1] = m_p.decode_all_mutation(test_stru, test_sec_1)
    p_mutant_mapper[test_sec_3] = m_p.decode_direct_mutation(test_stru, test_sec_3)
    mutants = m_p.combine_section_mutant(p_mutant_mapper)
    assert len(mutants) == 400
    for mut in mutants:
        assert Mutation(orig='LEU', target='ALA', chain_id='A', res_idx=10) in mut

def test_combine_section_mutant_many_to_many():
    """test the function works as expected"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_sec_1 = "a:[resi 253:all not self, resi 252:larger]"
    test_sec_2 = "r:2[resi 254 around 3:all not self]*5"
    # build per-section mutant mapper
    p_mutant_mapper = {}
    p_mutant_mapper[test_sec_1] = m_p.decode_all_mutation(test_stru, test_sec_1)
    p_mutant_mapper[test_sec_2] = m_p.decode_random_mutation(test_stru, test_sec_2)

    mutants = m_p.combine_section_mutant(p_mutant_mapper)
    assert len(mutants) == 2000
    for mut in mutants:
        assert len(mut) >= 2 and len(mut) <= 4
