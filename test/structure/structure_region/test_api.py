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
    """as name. TODO complete assert part"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "br. (resi 254 around 5)"
    )
    assert len(test_stru_region.atoms) == 335

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

def test_involved_residues_with_free_terminal():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 2+3+4")
    test_raw_region = stru_regi.StructureRegion(atoms=sele.atoms)
    result = test_raw_region.involved_residues_with_free_terminal()
    assert result["c_ter"] == [test_stru.get("A.4")]
    assert result["n_ter"] == [test_stru.get("A.2")]

def test_involved_residues_with_free_terminal_double():
    """test the case that contains residues that have both
    C and N ter exposed"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 2+5+7")
    test_raw_region = stru_regi.StructureRegion(atoms=sele.atoms)
    answer = set([
        test_stru.get("A.2"),
        test_stru.get("A.5"),
        test_stru.get("A.7"),
    ]
    )
    result = test_raw_region.involved_residues_with_free_terminal()
    assert set(result["c_ter"]) == answer
    assert set(result["n_ter"]) == answer

def test_involved_residues_with_free_terminal_ter():
    """test the case using chain terminal residues"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 1+253")
    test_raw_region = stru_regi.StructureRegion(atoms=sele.atoms)
    result = test_raw_region.involved_residues_with_free_terminal()
    assert result["c_ter"] == [test_stru.get("A.1")]
    assert result["n_ter"] == [test_stru.get("A.253")]

def test_get_net_charge():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "br. (resi 254 around 5) or resi 254"
    )
    assert test_stru_region.get_net_charge() == 0

def test_atoms_by_residue(): # TODO
    "as name. TODO finish this"
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin(
        {"H5J": (0,1)}
    )
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "br. (resi 254 around 5) or resi 254"
    )
    assert test_stru_region.atoms_by_residue
    assert False
