"""Testing the enzy_htp.structure.structure_region.capping.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-31
"""
import pytest
import os
from copy import deepcopy

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
import enzy_htp.core.file_system as fs
import enzy_htp.structure.structure_region.capping as capping
import enzy_htp.structure.structure_region as stru_regi
import enzy_htp.structure.structure_selection as stru_sele

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
WORK_DIR = f"{CURR_DIR}/../work_dir/"
sp = PDBParser()

def test_capping_with_residue_terminals(helpers):
    """as name. use an example region from selection.
    Just make sure not execptions is raised for now.
    TODO complete this test"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    sele = stru_sele.select_stru(
        test_stru, "resi 1+2+3+7+8")
    test_region = stru_regi.StructureRegion(atoms=sele.atoms)
    capping.capping_with_residue_terminals(test_region)
    # test
    atoms = test_region.atoms
    test_file = f"{WORK_DIR}test_capping_1.xyz"
    answer_file = f"{DATA_DIR}answer_capping_1.xyz"
    lines = [
        str(len(atoms)),
        "",]
    for at in atoms:
        lines.append(f"{at.element} {at.coord[0]} {at.coord[1]} {at.coord[2]}")
    fs.write_lines(test_file, lines)
    assert helpers.equiv_files(test_file, answer_file, consider_order=False)

    fs.safe_rm(test_file)

def test_residue_cap_deepcopy():
    """test the behavior of deepcopy for residue cap"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_res = test_stru.residues[0]
    test_res_2 = test_stru.residues[1]
    test_ch3 = capping.get_ch3_cap("cterm", test_res, test_res.atoms[0], test_res_2.atoms[0])

    ch3_copy = deepcopy(test_ch3)

    assert ch3_copy.link_residue is not test_res
    assert ch3_copy.link_residue.parent is None
    assert ch3_copy.plug_atom in ch3_copy.atoms
    assert ch3_copy.link_atom in ch3_copy.link_residue.atoms