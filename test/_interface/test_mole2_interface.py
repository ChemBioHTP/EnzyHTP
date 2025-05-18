"""Testing enzy_htp._interface.mole2_interface.py
Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Date: 2025-05-18
"""

from os import path
from typing import List
from enzy_htp import interface, _LOGGER, PDBParser
from enzy_htp._interface import Mole2Cavity

sp = PDBParser()
mole2 = interface.mole2

DATA_DIR = path.join(path.dirname(__file__), "data")
WORK_DIR = path.join(path.dirname(__file__), "work_dir")

def test_identify_cavities():
    """Test the `interface.mole2.identify_cavities` function."""
    pdb_filepath = path.join(DATA_DIR, "aclHMT-ETI-SAH_no-ETI.pdb")
    cavities: List[Mole2Cavity] = interface.mole2.identify_cavities(pdb_filepath, work_dir=WORK_DIR)
    
    has_target_cavity = False
    for i, cavity in enumerate(cavities):
        _LOGGER.info(f"Cavity {i+1}, Volume: {cavity.volume()}")
        if (cavity.volume() > 900):
            has_target_cavity = True
        else:
            pass
        continue
    assert has_target_cavity