"""Testing for the functions found in enzy_htp.preparation.protonate.py.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os

import pytest
import numpy as np
from enzy_htp import core as cc
from enzy_htp.core import file_system as fs
import enzy_htp.structure as struct
import enzy_htp.structure.structure_operation as stru_oper
from enzy_htp.preparation import protonate as prot


CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data/" 
WORK_DIR = f"{CURR_DIR}/work_dir/"
sp = struct.PDBParser()

# TODO(CJ): make testing utility files.
def equiv_files(fname1: str, fname2: str, width: int = None) -> bool:
    """Helper method to check if two files are exactly equivalent."""
    for l1, l2 in zip(fs.lines_from_file(fname1), fs.lines_from_file(fname2)):
        if width:
            l1 = l1[:width]
            l2 = l2[:width]

        if l1 != l2:
            print(f"'{l1}'")
            print(f"'{l2}'")
            return False
    return True

def test_protonate_stru():
    assert False

def test_protonate_peptide_with_pdb2pqr_no_metal():
    """test that protonate_peptide_with_pdb2pqr() works without exceptions"""
    test_pdb = f"{DATA_DIR}/1Q4T_cofactor_2chain_not_from_1.pdb"
    int_pdb = f"{WORK_DIR}/1Q4T_cofactor_2chain_not_from_1_in.pdb"
    int_pqr = f"{WORK_DIR}/1Q4T_cofactor_2chain_not_from_1_out.pdb"
    assert not os.path.exists(int_pdb)
    assert not os.path.exists(int_pqr)

    stru = sp.get_structure(test_pdb)
    prot.protonate_peptide_with_pdb2pqr(stru, 7.0, int_pdb, int_pqr)
    assert len(stru.atoms) == 4346
    assert len(stru.find_residue_name("HIS")) == 0
    assert len(stru.find_residue_name("HID")) == 10
    assert len(stru.find_residue_name("HIE")) == 2
    assert not os.path.exists(int_pqr)
    assert not os.path.exists(int_pdb)

def test_protonate_peptide_with_pdb2pqr_metal():
    """test that protonate_peptide_with_pdb2pqr() works without exceptions"""
    test_pdb = f"{DATA_DIR}/1NVG_metalcenter_noligand.pdb"
    int_pdb = f"{WORK_DIR}/1NVG_metalcenter_noligand_in.pdb"
    int_pqr = f"{WORK_DIR}/1NVG_metalcenter_noligand_out.pdb"
    # assert not os.path.exists(int_pdb)
    # assert not os.path.exists(int_pqr)

    stru = sp.get_structure(test_pdb)
    prot.protonate_peptide_with_pdb2pqr(stru, 7.0, int_pdb, int_pqr)
    assert list(map(lambda r:r.name ,stru.metalcenters[0].get_donor_mapper())) == ['GLU', 'CYM', 'CYM', 'CYM']
    assert len(stru.atoms) == 5353 # removed 5 atoms
    assert not os.path.exists(int_pqr)
    assert not os.path.exists(int_pdb)

def test_pdb2pqr_protonate_pdb_FAcD():
    """Making sure the protonate_pdb() method works for the FAcD enzyme system."""
    test_pdb = f"{DATA_DIR}/FAcD-FA-ASP_rmW.pdb"
    answer_pqr = f"{DATA_DIR}/FAcD-FA-ASP_rmW.pqr"
    work_pqr = f"{WORK_DIR}/FAcD-FA-ASP_rmW.pqr"
    assert not os.path.exists(work_pqr)
    prot.pdb2pqr_protonate_pdb(test_pdb, work_pqr)
    assert os.path.exists(work_pqr)
    assert equiv_files(answer_pqr, work_pqr)
    fs.safe_rm(work_pqr)
    assert not os.path.exists(work_pqr)

def test_pdb2pqr_protonate_pdb_4NKK():
    """Making sure the protonate_pdb() method works for the 4NKK enzyme system."""
    test_pdb = f"{DATA_DIR}/4NKK_clean.pdb"
    target_pqr = f"{DATA_DIR}/4NKK_clean.pqr"
    actual_pqr = f"{WORK_DIR}/4NKK_clean.pqr"
    assert not os.path.exists(actual_pqr)
    prot.pdb2pqr_protonate_pdb(test_pdb, actual_pqr)
    assert os.path.exists(actual_pqr)
    assert equiv_files(target_pqr, actual_pqr)
    fs.safe_rm(actual_pqr)
    assert not os.path.exists(actual_pqr)

def test_protonate_ligand_with_pybel():
    """test that protonate_ligand_with_pybel() works without exceptions"""
    test_pdb = f"{DATA_DIR}/FAcD-FA-ASP_rmW.pdb"

    stru = sp.get_structure(test_pdb)
    prot.protonate_ligand_with_pybel(stru, 7.0, int_pdb, int_pqr)
    assert len(stru.atoms) == 4346
    assert len(stru.find_residue_name("HIS")) == 0
    assert len(stru.find_residue_name("HID")) == 10
    assert len(stru.find_residue_name("HIE")) == 2

    fs.safe_rm(int_pqr)
    fs.safe_rm(int_pdb)
    assert not os.path.exists(int_pqr)
    assert not os.path.exists(int_pdb)
    

def test_protonate_ligand():
    """Ensuring that the protonate_ligand() method works."""
    # TODO(CJ): remove the old stuff we dont care about
    file_to_remove = list()
    ligand_dir = f"{WORK_DIR}/ligands/"
    base_pdb = f"{DATA_DIR}/FAcD-FA-ASP.pdb"
    structure: struct.Structure = struct.structure_from_pdb(base_pdb)
    ligand = structure.ligands[0]
    ligand.build(f"{WORK_DIR}/ligand.pdb")
    # assert not os.path.isdir( ligand_dir )
    protonated_ligand: struct.Ligand = prot.protonate_ligand(ligand, ligand_dir)
    assert protonated_ligand.net_charge == -1
    fs.safe_rmdir(ligand_dir)
    assert not os.path.isdir(ligand_dir)


def test_protonate_missing_elements_ligand():
    """Testing protonate_missing_elements() for PDB with a ligand included."""
    assert False


def test_protonate_missing_elements_metal():
    """Testing protonate_missing_elements() for PDB with a metal included."""
    assert False
