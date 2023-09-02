"""Testing the enzy_htp.structure.structure_operation.general.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-22
"""

import logging
import os

from enzy_htp.structure import PDBParser, Structure
import enzy_htp.structure.structure_operation as stru_oper
from enzy_htp.core import _LOGGER

_LOGGER.setLevel(logging.CRITICAL)

CURRDIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURRDIR}/../data/"
sp = PDBParser()


def test_remove_solvent():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)

    stru_oper.remove_solvent(stru)
    assert len(stru.solvents) == 0


def test_remove_empty_chain():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)

    stru[1].residues = []
    stru[2].residues = []  # be careful when delete element from a iterable
    stru_oper.remove_empty_chain(stru)
    assert len(stru) == 2


def test_remove_non_peptide():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    stru_oper.remove_non_peptide(stru)
    assert tuple(map(lambda x: x.name, stru)) == ("A", "B")


def test_update_residues():
    """test updating residues from a protonated structure"""
    pdb_file_path = f"{DATA_DIR}1Q4T_residue_update_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    new_pdb_file_path = f"{DATA_DIR}1Q4T_peptide_protonated.pdb"
    ref_stru: Structure = sp.get_structure(new_pdb_file_path)
    stru_oper.update_residues(stru, ref_stru)
    # this truncated the part with mem id
    assert str(stru)[37:] == ("Structure(\n"
                              "chains: (sorted, original [\'A\', \'B\', \'C\', \'D\'])\n"
                              "    A(peptide): residue: 10-151 atom_count: 2127\n"
                              "    B(peptide): residue: 12-151 atom_count: 2103\n"
                              "    C(ligand): residue: 370-370 atom_count: 58\n"
                              "    D(ligand): residue: 371-371 atom_count: 58\n"
                              ")")


def test_align_atom_order_in_each_residue():
    """test updating the atom order in each residue of an enzyme"""
    pdb_file_path = f"{DATA_DIR}KE_07_R7_2_S_mut.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    pdb_file_path_2 = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    stru_2: Structure = sp.get_structure(pdb_file_path_2)
    stru_oper.align_atom_order_in_each_residue(stru, stru_2)
    
    for new_res, old_res in zip(stru.residues, stru_2.residues):
        if new_res.idx == 154:
            assert new_res.name == "TRP"
            assert len(new_res.atoms) == 24
        else:
            for new_atom, old_atom in zip(new_res.atoms, old_res.atoms):
                assert new_atom.name == old_atom.name and new_atom.idx == old_atom.idx

    