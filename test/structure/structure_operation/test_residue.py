"""Testing the enzy_htp.structure.structure_operation.residue.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-03-20
"""

import logging
import os

import pytest
from enzy_htp.structure import PDBParser, Structure
import enzy_htp.structure.structure_operation as stru_oper
from enzy_htp.core import _LOGGER

_LOGGER.setLevel(logging.CRITICAL)

CURRDIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURRDIR}/../data/"
sp = PDBParser()


def test_deprotonate_residue():
    """test deprotonate_residue considered CYS on SG and HID on NE2 and ND1"""
    pdb_file_path = f"{DATA_DIR}1Q4T_peptide_protonated.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    # CYS A69
    residue = stru["A"].find_residue_idx(69)
    target_atom = residue.find_atom_name("SG")
    stru_oper.deprotonate_residue(residue, target_atom)
    assert residue.name == "CYM"
    assert len(residue) == 10
    # HID A21 (already no H)
    ## case where no need to deproton
    residue = stru["A"].find_residue_idx(21)
    target_atom = residue.find_atom_name("NE2")
    stru_oper.deprotonate_residue(residue, target_atom)
    assert residue.name == "HID"
    assert len(residue) == 17


@pytest.mark.TODO  # need add a method to complete atoms
def test_deprotonate_residue_switch():
    """test deprotonate_residue considered CYS on SG and HID on NE2 and ND1"""
    pdb_file_path = f"{DATA_DIR}1Q4T_peptide_protonated.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    # manually select the case where it is a switch
    residue = stru["A"].find_residue_idx(21)  # HID 21
    target_atom = residue.find_atom_name("ND1")
    stru_oper.deprotonate_residue(residue, target_atom)

    assert residue.name == "HIE"
    assert len(residue) == 17
    # TODO add test with completion of LYS treatment etc.


def test_remove_side_chain_atom():
    """test function works as expected"""
    pdb_file_path = f"{DATA_DIR}1Q4T_peptide_protonated.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    # manually select the case where it is a switch
    test_residue = stru.residues[5]
    stru_oper.remove_side_chain_mutating_atom(test_residue, "ASP")
    assert set(test_residue.atom_name_list) == set(["N", "CA", "C", "O", "CB", "H", "HA"])
    test_residue = stru.residues[7]
    stru_oper.remove_side_chain_mutating_atom(test_residue, "GLY")
    assert set(test_residue.atom_name_list) == set(["N", "CA", "C", "O", "H"])
