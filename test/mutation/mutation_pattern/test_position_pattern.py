"""Testing the enzy_htp.mutation.mutation_pattern.position_pattern.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-17
"""

import os

from enzy_htp import PDBParser
import enzy_htp.mutation.mutation_pattern.general as m_p

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

def test_decode_position_pattern():
    """test the function use a made up position_pattern for KE"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 254 around 4"

    assert m_p.decode_position_pattern(test_stru, test_pattern) == [
        ('A', 224), ('A', 202), ('A', 9), ('A', 201), ('A', 48), ('A', 11), ('A', 169),
        ('A', 222), ('A', 101), ('A', 128), ('A', 50)]

def test_decode_position_pattern_dimer():
    """test the function use a made up position_pattern for puo"""
    test_pdb = f"{DATA_DIR}puo_put.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 902 around 4"

    assert m_p.decode_position_pattern(test_stru, test_pattern) == [
        ('A', 205), ('A', 394), ('A', 431), ('A', 59), ('A', 323), ('A', 171),
        ('A', 340), ('A', 172), ('C', 901), ('A', 325)]

