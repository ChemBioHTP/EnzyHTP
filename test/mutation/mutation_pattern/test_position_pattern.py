"""Testing the enzy_htp.mutation.mutation_pattern.position_pattern.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-17
"""

import os

from enzy_htp import PDBParser
import enzy_htp.mutation.mutation_pattern.api as m_p

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()


def test_decode_position_pattern():
    """test the function use a made up position_pattern for KE"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 254 around 4"

    assert set(m_p.decode_position_pattern(test_stru, test_pattern)) == set([
        ('A', 224), ('A', 202), ('A', 9), ('A', 201), ('A', 48), ('A', 11), ('A', 169),
        ('A', 222), ('A', 101), ('A', 128), ('A', 50)
    ])


def test_decode_position_pattern_dimer():
    """test the function use a made up position_pattern for puo"""
    test_pdb = f"{DATA_DIR}puo_put.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 267 around 5"
    assert set(m_p.decode_position_pattern(test_stru, test_pattern)) == set([
        ('A', 273), ('A', 268), ('A', 263), ('A', 389), ('A', 290), ('A', 269),
        ('A', 264), ('A', 390), ('A', 270), ('A', 265), ('A', 391), ('A', 271),
        ('A', 266), ('A', 388), ('B', 719)
    ]) # note ligand wont be included
