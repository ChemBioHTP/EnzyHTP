"""Testing the enzy_htp.structure.structure_constraint.api.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-25
"""
import pytest
import os

from enzy_htp import PDBParser
import enzy_htp.structure.structure_constraint as stru_cons

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

def test_create_backbone_freeze():
    "test function works as expected"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_bb_freeze = stru_cons.create_backbone_freeze(test_stru)
    assert test_bb_freeze.params["amber"] == {"restraint_wt": 2.0}
    assert test_bb_freeze.constraint_type == "backbone_freeze"
    assert test_bb_freeze.atom_names == {"C", "CA", "N"}
