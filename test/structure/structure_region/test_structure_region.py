"""Testing the enzy_htp.structure.structure_region.structure_region.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-02
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
    test_stru.assign_ncaa_chargespin(
        {"H5J": (0,1)}
    )
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "br. (resi 254 around 5) or resi 254"
    )
    init_charge(test_stru_region)
    assert test_stru_region.get_net_charge() == 0

def test_get_spin():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin(
        {"H5J": (0,1)}
    )
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "br. (resi 254 around 5) or resi 254"
    )
    # init_spin(test_stru_region)
    assert test_stru_region.get_spin() == 1

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

def test_atoms_by_residue_simple(): # TODO
    "as name. TODO finish this"
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin(
        {"H5J": (0,1)}
    )
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "resi 2+3+4"
    )
    # print(*map(lambda x: list(x), test_stru_region.atoms_by_residue.values()), sep="\n")

def test_atoms_from_geom(helpers):
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "resi 2+3+4", nterm_cap="CH3", cterm_cap="CH3"
    )
    test_geom = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S_geom_1.pdb")

    result = test_stru_region.atoms_from_geom(test_geom)

    # test TODO refine this after xyz interface
    atoms = result
    test_file = f"{WORK_DIR}test_atom_from_geom_1.xyz"
    answer_file = f"{DATA_DIR}answer_atom_from_geom_1.xyz"
    lines = [
        str(len(atoms)),
        "",]
    for at in atoms:
        lines.append(f"{at.element} {at.coord[0]} {at.coord[1]} {at.coord[2]}")
    fs.write_lines(test_file, lines)
    assert helpers.equiv_files(test_file, answer_file, consider_order=False)

    fs.safe_rm(test_file)

def test_topology():
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "all"
    )
    test_stru_2 = test_stru_region.topology
    assert test_stru.is_same_topology(test_stru_2)

def test_convert_to_structure():
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J": (0,1)})
    test_stru_region = stru_regi.create_region_from_selection_pattern(
        test_stru, "all"
    )
    test_stru_2 = test_stru_region.convert_to_structure()
    assert test_stru.is_same_topology(test_stru_2)

    test_stru_3 = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_3.assign_ncaa_chargespin({"H5J": (0,1)})
    test_stru_region_2 = stru_regi.create_region_from_full_stru(
        test_stru_3
    )
    test_stru_4 = test_stru_region_2.convert_to_structure()
    assert test_stru_3.is_same_topology(test_stru_4)

