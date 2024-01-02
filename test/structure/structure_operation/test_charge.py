"""Testing the functionality of the enzy_htp.structure.charge.py
Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-01-01
"""
import os
import pytest
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Chain, Residue, PDBParser, Structure
import enzy_htp.structure.structure_operation.charge as chrg

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"
sp = PDBParser()

def test_init_charge_caa():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_res = test_stru.residues[0]

    chrg.init_charge(test_res)

