"""Testing the enzy_htp.structure.structure_region.residue_caps.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-31
"""
import pytest
import os
from copy import deepcopy
from typing import List

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
import enzy_htp.core.file_system as fs
import enzy_htp.structure.structure_region.capping as capping
import enzy_htp.structure.structure_region as stru_regi
import enzy_htp.structure.structure_selection as stru_sele
from enzy_htp.structure.structure_region.residue_caps import SUPPORTED_CAPS, ResidueCap, CH3Cap 

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
WORK_DIR = f"{CURR_DIR}/../work_dir/"
sp = PDBParser()

def test_residue_cap_deepcopy():
    """test the behavior of deepcopy for residue cap"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_res = test_stru.residues[0]
    test_res_2 = test_stru.residues[1]
    test_ch3 = CH3Cap(
        test_res,
        test_res.find_atom_name('C'),
        test_res_2.find_atom_name('N'),
        'cterm' 
    )
    ch3_copy = deepcopy(test_ch3)

    assert ch3_copy.link_residue is not test_res
    assert ch3_copy.link_residue.parent is None
    assert ch3_copy.link_atom in ch3_copy.link_residue.atoms
