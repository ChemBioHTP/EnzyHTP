"""Testing the enzy_htp.structure.structure_selection.general.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-15
"""

import pytest
import os

from enzy_htp import PDBParser
from enzy_htp.structure.structure_selection import select_stru

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

@pytest.mark.interface
def test_select_stru():
    "test function works as expected"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 289 around 4"
    select_stru(test_stru, test_pattern)
