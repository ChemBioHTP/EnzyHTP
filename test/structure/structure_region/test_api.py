"""Testing the enzy_htp.structure.structure_region.api.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-31
"""
import pytest
import os

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
import enzy_htp.structure.structure_region as stru_regi
import enzy_htp.structure.structure_selection as stru_sele

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

def test_create_region_from_selection_pattern():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        "br. (resi 254 around 5)", test_stru
    )

def test_is_whole_residue_only():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    sele = stru_sele.select_stru(
        test_stru, "br. (resi 254 around 5)")
    test_raw_region = stru_regi.StructureRegion(atoms=sele.atoms)
    assert test_raw_region.is_whole_residue_only()

def test_is_whole_residue_only_false():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 254 around 5")
    test_raw_region = stru_regi.StructureRegion(atoms=sele.atoms)
    assert not test_raw_region.is_whole_residue_only()
