"""Testing the Structure() class. 
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""
import os
import pytest
import numpy as np
from copy import deepcopy
import enzy_htp.chemical as chem
from enzy_htp.core import file_system as fs
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp.structure import (
    Structure,
    Residue,
    Chain,
    Atom,
)

#TODO(CJ): add tests for Structure.build_ligands()
CURR_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.dirname(CURR_DIR)
DATA_DIR = f"{CURR_DIR}/data/"
sp = PDBParser()


def test_structure_ctor_bad_input():
    """Testing that the Structure() ctor fails when given duplicate chains."""
    Atom_1 = Atom({"atom_name": "X", "x_coord": 1, "y_coord": 2, "z_coord": 3})
    Atom_2 = Atom({"atom_name": "X", "x_coord": 3, "y_coord": 2, "z_coord": 1})
    chain1 = Chain("A", [Residue(1, "XXX", [Atom_1])])
    chain2 = Chain("B", [Residue(2, "XXX", [Atom_2])])
    chain3 = Chain("A", [Residue(3, "XXX", [Atom_1])])

    with pytest.raises(SystemExit) as exe:
        struct = Structure([chain1, chain2, chain3])

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1


def test_residues():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)

    assert len(stru.residues) == 580


def test_ligands():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert len(stru.ligands) == 2
    assert stru.ligands[0].name == "4CO"


def test_solvents():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert len(stru.solvents) == 296
    assert stru.solvents[0].name == "HOH"


def test_peptides():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert len(stru.peptides) == 2
    assert tuple(map(lambda x: x.name, stru.peptides)) == ("A", "B")


def test_sequence():
    pdb_file_path = f"{DATA_DIR}1Q4T_ligand_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    target_seq: str = "ATGGNLPDVASHYPVAYEQTLDGTVGFVIDEMTPERATASVEVTDTLRQRWGLVHGGAYCALAEMLATEATVAVVHEKGMMAVGQSNHTSFFRPVKEGHVRAEAVRIHAGSTTWFWDVSLRDDAGRLCAVSSMSIAVRPRRD"
    assert stru.sequence["A"] == target_seq


def test_is_idx_subset_subset():
    """test the checker for idx subset with a handmade subset stru"""
    pdb_file_path = f"{DATA_DIR}12E8_small_four_chain.pdb"
    subset_pdb_file_path = f"{DATA_DIR}12E8_small_four_chain_subset.pdb"

    stru: Structure = sp.get_structure(pdb_file_path)
    target_stru: Structure = sp.get_structure(subset_pdb_file_path)
    assert stru.is_idx_subset(target_stru)


def test_is_idx_subset_non_res_subset():
    """test the checker for idx subset with a handmade non-subset (residue) stru"""
    pdb_file_path = f"{DATA_DIR}12E8_small_four_chain.pdb"
    non_subset_pdb_file_path = f"{DATA_DIR}12E8_small_four_chain_non_seq_subset.pdb"

    stru: Structure = sp.get_structure(pdb_file_path)
    target_stru: Structure = sp.get_structure(non_subset_pdb_file_path)
    assert not stru.is_idx_subset(target_stru)


def test_is_idx_subset_non_ch_subset():
    """test the checker for idx subset with a handmade non-subset (chain id) stru"""
    pdb_file_path = f"{DATA_DIR}12E8_small_four_chain.pdb"
    non_subset_pdb_file_path = f"{DATA_DIR}12E8_small_four_chain_non_chain_subset.pdb"

    stru: Structure = sp.get_structure(pdb_file_path)
    target_stru: Structure = sp.get_structure(non_subset_pdb_file_path)
    assert not stru.is_idx_subset(target_stru)


def equiv_files(fname1: str, fname2: str, width: int = None) -> bool:
    """Helper method to check if two files are exactly equivalent."""

    def find_first_diff(l1, l2):
        for idx, (ch1, ch2) in enumerate(zip(l1, l2)):
            if ch1 != ch2:
                return idx + 1
        return -1

    for line_idx, (l1, l2) in enumerate(
            zip(fs.lines_from_file(fname1), fs.lines_from_file(fname2))):
        if width:
            l1 = l1[:width]
            l2 = l2[:width]
            if len(l2) < width:
                l2 += " " * (width - len(l2))
            if len(l1) < width:
                l1 += " " * (width - len(l1))

        first_diff = find_first_diff(l1, l2)
        if first_diff != -1:
            diff = [" "] * min(len(l1), len(l2))
            diff[first_diff] = "^"

            print(f"'{l1}'")
            print(f"'{l2}'")
            print("".join(diff))
            print(
                f"Difference encountered on line {line_idx}. '{fname1}' and '{fname2}' are NOT equivalent"
            )
            return False
    return True


def test_atoms():  # TODO(shaoqz) wait for test
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    struct: Structure = PDBParser().get_structure(TEST_FILE)
    assert struct.atoms
