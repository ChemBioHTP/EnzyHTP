"""Testing the enzy_htp.structure.structure_region.api.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-31
"""
import pytest
import os

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
import enzy_htp.core.file_system as fs
import enzy_htp.structure.structure_region as stru_regi
import enzy_htp.structure.structure_selection as stru_sele
from enzy_htp.structure.structure_enchantment import init_charge

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
WORK_DIR = f"{CURR_DIR}/../work_dir/"
sp = PDBParser()

def test_create_region_from_selection_pattern():
    """as name. TODO complete assert part"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "br. (resi 254 around 5)", nterm_cap="CH3", cterm_cap="CH3"
    )
    assert len(test_stru_region.atoms) == 335

