"""Testing for the functions found in enzy_htp.preparation.protonate.py.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os

import pytest
import numpy as np
from enzy_htp import core as cc
from enzy_htp.core import file_system as fs
from enzy_htp import config
from enzy_htp.core.logger import _LOGGER
import logging
import enzy_htp.structure as struct
import enzy_htp.structure.structure_operation as stru_oper
from enzy_htp.preparation import protonate as prot

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data/"
WORK_DIR = f"{CURR_DIR}/work_dir/"
config["system.SCRATCH_DIR"] = WORK_DIR
sp = struct.PDBParser()


# TODO(CJ): make testing utility files.
def equiv_files(fname1: str, fname2: str, width: int = None, skip_frist: bool = False) -> bool:
    """Helper method to check if two files are exactly equivalent."""
    first = 1
    for l1, l2 in zip(fs.lines_from_file(fname1), fs.lines_from_file(fname2)):
        if skip_frist and first:
            first = 0
            continue
        if width:
            l1 = l1[:width]
            l2 = l2[:width]

        if l1 != l2:
            print(f"'{l1}'")
            print(f"'{l2}'")
            return False
    return True


ligand_answer_list_1Q4T = [
    "N1A", "C2A", "N3A", "C4A", "C5A", "C6A", "N6A", "N7A", "C8A", "N9A", "C1D", "C2D", "O2D", "C3D", "O3D", "P3D", "O7A", "O8A", "O9A",
    "C4D", "O4D", "C5D", "O5D", "P1A", "O1A", "O2A", "O3A", "P2A", "O4A", "O5A", "O6A", "CBP", "CCP", "CDP", "CEP", "CAP", "OAP", "C9P",
    "O9P", "N8P", "C7P", "C6P", "C5P", "O5P", "N4P", "C3P", "C2P", "S1P", "O1B", "C1B", "C2B", "C3B", "C4B", "C5B", "O2B", "C6B", "C7B",
    "CB", "H", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15", "H16", "H17", "H18", "H19",
    "H20", "H21", "H22", "H23", "H24", "H25", "H26", "H27", "H28", "H29", "H30", "H31", "H32", "H33", "H34", "H35", "H36", "H37"
]


def test_protonate_stru():
    """test the interface function mimicing how it will be used.
    Use 1Q4T as model structure. Use default settins"""
    test_pdb = f"{DATA_DIR}/1Q4T_cofactor_2chain_not_from_1.pdb"

    stru = sp.get_structure(test_pdb)
    prot.protonate_stru(stru, 7.0, protonate_ligand=True)
    #peptide
    assert len(stru.atoms) == 4422
    assert len(stru.find_residue_name("HIS")) == 0
    assert len(stru.find_residue_name("HID")) == 10
    assert len(stru.find_residue_name("HIE")) == 2
    #ligand
    assert stru.ligands[0].atom_name_list == ligand_answer_list_1Q4T
    assert stru.ligands[0].idx == 370
    assert stru.ligands[0].name == "4CO"
    assert stru.ligands[1].atom_name_list == ligand_answer_list_1Q4T
    assert stru.ligands[1].idx == 371
    assert stru.ligands[1].name == "4CO"


def test_protonate_stru_metal():
    """test the interface function mimicing how it will be used.
    Use 1NVG as model structure. Use default settings"""
    test_pdb = f"{DATA_DIR}/1NVG_metalcenter_noligand.pdb"

    stru = sp.get_structure(test_pdb)
    prot.protonate_stru(stru, 7.0, protonate_ligand=True)
    #peptide
    assert len(stru.atoms) == 5353
    assert len(stru.find_residue_name("HIS")) == 0
    assert len(stru.find_residue_name("HID")) == 6
    assert len(stru.find_residue_name("HIE")) == 1
    #metal
    assert len(stru.find_residue_name("CYM")) == 5


def test_protonate_peptide_with_pdb2pqr_no_metal():
    """test that protonate_peptide_with_pdb2pqr() works without exceptions"""
    test_pdb = f"{DATA_DIR}/1Q4T_cofactor_2chain_not_from_1.pdb"
    int_pdb = f"{WORK_DIR}/1Q4T_cofactor_2chain_not_from_1_in.pdb"
    int_pqr = f"{WORK_DIR}/1Q4T_cofactor_2chain_not_from_1_out.pdb"
    fs.safe_mkdir(WORK_DIR)

    fs.safe_rm(int_pdb)
    fs.safe_rm(int_pqr)

    assert not os.path.exists(int_pdb)
    assert not os.path.exists(int_pqr)

    stru = sp.get_structure(test_pdb)
    prot.protonate_peptide_with_pdb2pqr(stru, 7.0, int_pdb, int_pqr)
    assert len(stru.atoms) == 4346
    assert len(stru.find_residue_name("HIS")) == 0
    assert len(stru.find_residue_name("HID")) == 10
    assert len(stru.find_residue_name("HIE")) == 2
    assert os.path.exists(int_pqr)  # if specified, it should be preserved
    assert os.path.exists(int_pdb)
    fs.safe_rm(int_pqr)
    fs.safe_rm(int_pdb)


def test_protonate_peptide_with_pdb2pqr_remove_temp():
    """test that protonate_peptide_with_pdb2pqr() works removing temp file"""
    test_pdb = f"{DATA_DIR}/1Q4T_cofactor_2chain_not_from_1.pdb"
    scratch_dir = config["system.SCRATCH_DIR"]
    int_pdb = f"{scratch_dir}/protonate_peptide_with_pdb2pqr_input.pdb"
    int_pqr = f"{scratch_dir}/protonate_peptide_with_pdb2pqr_output.pdb"

    fs.safe_rm(int_pdb)
    fs.safe_rm(int_pqr)

    assert not os.path.exists(int_pdb)
    assert not os.path.exists(int_pqr)

    stru = sp.get_structure(test_pdb)
    prot.protonate_peptide_with_pdb2pqr(stru, 7.0)
    if _LOGGER.level > 10:
        assert not os.path.exists(int_pqr)
        assert not os.path.exists(int_pdb)


def test_protonate_peptide_with_pdb2pqr_metal():
    """test that protonate_peptide_with_pdb2pqr() works without exceptions"""
    test_pdb = f"{DATA_DIR}/1NVG_metalcenter_noligand.pdb"
    int_pdb = f"{WORK_DIR}/1NVG_metalcenter_noligand_in.pdb"
    int_pqr = f"{WORK_DIR}/1NVG_metalcenter_noligand_out.pdb"
    fs.safe_mkdir(WORK_DIR)
    # assert not os.path.exists(int_pdb)
    # assert not os.path.exists(int_pqr)

    stru = sp.get_structure(test_pdb)
    prot.protonate_peptide_with_pdb2pqr(stru, 7.0, int_pdb, int_pqr)
    assert list(map(lambda r: r.name, stru.metalcenters[0].get_donor_mapper())) == ["GLU", "CYM", "CYM", "CYM"]
    assert len(stru.atoms) == 5353  # removed 5 atoms
    assert os.path.exists(int_pqr)  # if specified, it should be preserved
    assert os.path.exists(int_pdb)
    fs.safe_rm(int_pqr)
    fs.safe_rm(int_pdb)


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
    """Making sure the protonate_pdb() method works for the 4NKK enzyme system.
    Note: there could be errors of this related to different PDB2PQR version"""
    test_pdb = f"{DATA_DIR}/4NKK_clean.pdb"
    target_pqr = f"{DATA_DIR}/4NKK_clean.pqr"
    target_pqr_old = f"{DATA_DIR}/4NKK_clean_old.pqr"
    actual_pqr = f"{WORK_DIR}/4NKK_clean.pqr"
    assert not os.path.exists(actual_pqr)
    prot.pdb2pqr_protonate_pdb(test_pdb, actual_pqr)
    assert os.path.exists(actual_pqr)
    assert equiv_files(target_pqr, actual_pqr) or equiv_files(target_pqr_old, actual_pqr)
    fs.safe_rm(actual_pqr)
    assert not os.path.exists(actual_pqr)


def test_protonate_ligand_with_pybel():
    """test that protonate_ligand_with_pybel() works without exceptions
    verify the protonation with manual curated atom list  and residue name/idx.
    also test normal behavior of removing int_dir"""
    test_pdb = f"{DATA_DIR}/FAcD-FA-ASP_rmH.pdb"
    int_lig_dir = f"{WORK_DIR}/ligand_test"
    assert not os.path.isdir(int_lig_dir)

    stru = sp.get_structure(test_pdb)

    prot.protonate_ligand_with_pybel(stru, 7.0, int_lig_dir)

    assert stru.ligands[0].atom_name_list == ["C", "F", "O", "CH3", "OXT", "H", "H1"]
    assert stru.ligands[0].idx == 298
    assert stru.ligands[0].name == "FAH"

    assert not os.path.isdir(int_lig_dir)

def test_protonate_ligand_with_pybel_7si9():
    """Test the protonation to ligands in a structure instance"""

    file_prefix = '7SI9_E166Q_protonated'
    pdb_file = f"{DATA_DIR}{file_prefix}.pdb"
    int_ligand_file_dir=f"{WORK_DIR}"
    
    stru = sp.get_structure(pdb_file)
    prot.protonate_ligand_with_pybel(stru=stru, int_ligand_file_dir=int_ligand_file_dir)
    
    for atom in stru.ligands[0].atoms:
        element_from_name = atom.name[:len(atom.element)]   # The name of the atom should start with the element symbol.
        assert element_from_name.upper() == atom.element.upper()
    

def test_protonate_ligand_with_pybel_2_ligand():
    """test that protonate_ligand_with_pybel() works without exceptions using 1Q4T
    which contains 2 complex ligands verify the protonation with manual curated atom
    list and residue name/idx. also test normal behavior removing int_dir"""
    test_pdb = f"{DATA_DIR}/1Q4T_cofactor_2chain_not_from_1.pdb"
    int_lig_dir = f"{WORK_DIR}/ligand_test"

    assert not os.path.isdir(int_lig_dir)

    stru = sp.get_structure(test_pdb)
    prot.protonate_ligand_with_pybel(stru, 7.0, int_lig_dir)

    assert stru.ligands[0].atom_name_list == ligand_answer_list_1Q4T
    assert stru.ligands[0].idx == 370
    assert stru.ligands[0].name == "4CO"
    assert stru.ligands[1].atom_name_list == ligand_answer_list_1Q4T
    assert stru.ligands[1].idx == 371
    assert stru.ligands[1].name == "4CO"

    assert not os.path.isdir(int_lig_dir)


def test_pybel_protonate_pdb_ligand():
    """test this function run without any error. use a manually verified file as answer"""
    ligand_path = f"{DATA_DIR}/ligand_test_FAH.pdb"
    out_ligand_path = f"{WORK_DIR}/ligand_test_FAH_pybel.pdb"
    answer_ligand_path = f"{DATA_DIR}/ligand_test_FAH_pybel.pdb"
    assert not os.path.isdir(out_ligand_path)

    prot.pybel_protonate_pdb_ligand(ligand_path, out_ligand_path)

    assert equiv_files(out_ligand_path, answer_ligand_path, skip_frist=True)
    fs.safe_rm(out_ligand_path)
    assert not os.path.isdir(out_ligand_path)


def test_pybel_protonate_pdb_ligand_4CO():
    """test this function run without any error. use a manually verified file as answer"""
    ligand_path = f"{DATA_DIR}/ligand_test_4CO.pdb"
    out_ligand_path = f"{WORK_DIR}/ligand_test_4CO_pybel.pdb"
    answer_ligand_path = f"{DATA_DIR}/ligand_test_4CO_pybel.pdb"
    assert not os.path.isdir(out_ligand_path)

    prot.pybel_protonate_pdb_ligand(ligand_path, out_ligand_path)

    # assert equiv_files(out_ligand_path, answer_ligand_path, skip_frist=True)
    # fs.safe_rm(out_ligand_path)
    # assert not os.path.isdir(out_ligand_path)

def test_fix_pybel_output():
    """test an exception case for fixing pybel output"""
    ligand_path = f"{DATA_DIR}ligand_test_HEZ_pybel_badname.pdb"
    ref_ligand_path = f"{DATA_DIR}ligand_test_HEZ.pdb"
    out_ligand_path = f"{WORK_DIR}ligand_test_HEZ_pybel.pdb"
    answer_ligand_path = f"{DATA_DIR}ligand_test_HEZ_pybel.pdb"
    assert not os.path.isdir(out_ligand_path)

    prot._fix_pybel_output(ligand_path, out_ligand_path, ref_ligand_path)

    assert equiv_files(out_ligand_path, answer_ligand_path, skip_frist=True)
    fs.safe_rm(out_ligand_path)
    assert not os.path.isdir(out_ligand_path)

def test_fix_pybel_output_4WI(caplog):
    """Test protonated ligand (pybel) name fixing."""
    file_prefix = 'ligand_test_4WI_protonated'
    pdb_file = f"{DATA_DIR}{file_prefix}.pdb"
    intermediate_pdb_file = f"{DATA_DIR}{file_prefix}_pybel_badname.pdb"
    out_pdb_file = f"{WORK_DIR}{file_prefix}_pybel.pdb"

    with pytest.raises(ValueError) as e:    # The value error is expected to be raised when unexpected Hydrogen atom(s) is detected.
        _LOGGER.propagate = True
        prot._fix_pybel_output(pdb_path=intermediate_pdb_file, out_path=out_pdb_file, ref_name_path=pdb_file)
    assert 'ref_name_path' in caplog.text
    # TODO (Zhong): The best way to do this in the _fix_pybel_output is: 
    # detect if there are any hydrogen interspersed in the middle of the file,
    # if the hydrogen atoms are all at the end of the file, it will not cause a problem.

def test_fix_pybel_output_4WI_H_end():
    """Test protonated ligand (pybel) name fixing.
    All the hydrogen atoms in the file of `ref_name` is at the end.
    The output file should be generated, pass the assert check and then be deleted.
    """
    file_prefix = 'ligand_test_4WI_protonated'
    pdb_file = f"{DATA_DIR}{file_prefix}_H_end.pdb"
    intermediate_pdb_file = f"{DATA_DIR}{file_prefix}_pybel_badname.pdb"
    out_pdb_file = f"{WORK_DIR}{file_prefix}_pybel.pdb"
    prot._fix_pybel_output(pdb_path=intermediate_pdb_file, out_path=out_pdb_file, ref_name_path=pdb_file)
    assert os.path.isfile(out_pdb_file)
    fs.clean_temp_file_n_dir([out_pdb_file])


def test_protonate_stru_imputed():
    """Testing the protonate_stru() method for a structure that has been imputed.
    TODO: those files are missing @CJ"""
    pdb_file = f"{DATA_DIR}/5k7u_imputed.pdb"
    stru = sp.get_structure(pdb_file)
    prot.protonate_stru(stru)

    assert stru

    pdb_file = f"{DATA_DIR}/6lfe_imputed.pdb"
    stru = sp.get_structure(pdb_file)
    prot.protonate_stru(stru)

    assert stru
