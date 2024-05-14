"""Testing the enzy_htp.structure.Residue() class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import copy
import itertools
import os
import pytest
import pandas as pd
from typing import List
from copy import deepcopy
from collections import defaultdict
from biopandas.pdb import PandasPdb
from enzy_htp.core.exception import ResidueDontHaveAtom

from enzy_htp.structure import Atom, Residue, Chain
from enzy_htp.structure.structure_io import PDBParser
from enzy_htp.structure.structure_enchantment import init_connectivity
import enzy_htp.structure.structure_operation as so

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"

TEST_PDB_FILE = f"{CURR_DIR}/data/3NIR.pdb"
RESIDUES = so.remove_solvent(PDBParser().get_structure(TEST_PDB_FILE)).residues

def test_has_hydrogens_FAH():
    """Test if hydrogen atoms can be properly detected."""
    file_prefix = 'FAcD-FA-ASP_rmW'
    pdb_file = f"{DATA_DIR}{file_prefix}.pdb"
    stru = PDBParser().get_structure(pdb_file)
    ligand = stru.ligands[0]
    assert ligand.has_hydrogens()

def test_has_hydrogens_4CO():
    """Test whether a ligand that does not contain a hydrogen atom can get the correct return value (False is expected)."""
    file_prefix = '1Q4T_residue_update_test'
    pdb_file = f"{DATA_DIR}{file_prefix}.pdb"
    stru = PDBParser().get_structure(pdb_file)
    ligand = stru.ligands[0]
    assert not ligand.has_hydrogens()
    

def test_deepcopy():
    """test the hehavior of copy.deepcopy on Residue() under a Structure()
    context"""
    stru = PDBParser().get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    res_list = stru[0][2:6] + stru[1][2:6]  # target for deepcopy
    new_list = deepcopy(res_list)
    # ensure the list is new
    assert id(res_list) != id(new_list)
    # ensure every residue is new
    for i, j in zip(res_list, new_list):
        assert id(i) != id(j)
    # ensure every parents is None
    assert all(i.parent is None for i in new_list)
    # ensure every atom in the stru is new and pointing to correct residue
    old_atom_list = itertools.chain.from_iterable(i.atoms for i in res_list)
    new_atom_list = itertools.chain.from_iterable(i.atoms for i in new_list)
    for i, j in zip(old_atom_list, new_atom_list):
        assert i.idx == j.idx
        assert id(i) != id(j)
    # ensure parent is containing the same new children
    for i in new_list:
        assert all(j.parent is i for j in i)


def test_find_atom_name():
    stru = PDBParser().get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    residue = stru["H"][0]
    assert residue.find_atom_name("N").name == "N"
    with pytest.raises(ResidueDontHaveAtom) as exe:
        residue.find_atom_name("X")
        assert exe.residue is residue
        assert exe.atom_name == "X"


#TODO recover tests
#TODO(CJ): add tests for the name getter
@pytest.mark.TODO
def make_residues(pdbname: str) -> List[Residue]:
    """Helper method that retrieves a list of residues from a PDB file."""
    holder = defaultdict(list)
    p_df: pd.DataFrame = PandasPdb().read_pdb(pdbname).df["ATOM"]

    for i, row in p_df.iterrows():
        atom: Atom = Atom.from_biopandas(**row)
        holder[atom.residue_key()].append(atom)

    residues = []
    for rkey, atoms in holder.items():
        residues.append(Residue(rkey, atoms))
    return residues


def test_check_all_canonical():
    """Checking that all the loaded Residue()'s are canonical amino acids and not metals, solvents of rd_non_ligands"""
    for rr in RESIDUES:
        print(rr)
        #assert rr.is_canonical()
        #assert not rr.is_metal()
        #assert not rr.is_solvent()
    assert not rr.is_trash()


def test_sequence_equivalent_return_true():
    """Testing that the Residue.is_sequence_equivalent() method returns True when desired."""
    r1 = Residue(1, "A", list(), Chain('A', list()))
    r2 = Residue(1, "A", list(), Chain('A', list()))
    assert r1.is_sequence_eq(r2)
    r3 = Residue(1, "A", list(), Chain('A', list()))
    r3.name = "B"
    assert r1.is_sequence_eq(r3)
    r3.chain.name = "A"
    assert r1.is_sequence_eq(r3)
    r3.chain.name = "B"
    assert not r1.is_sequence_eq(r3)


def test_sequence_equivalent_return_false():
    """Testing that the Residue.is_sequence_equivalent() method returns False when desired."""
    r1 = Residue(1, "A", list(), Chain('A', list()))
    r2 = Residue(3, "B", list(), Chain('A', list()))
    assert not r1.is_sequence_eq(r2)
    r3 = Residue(1, "C", list(), Chain('B', list()))
    assert not r1.is_sequence_eq(r3)


# TODO(CJ) add in some tests for other types of residues


def test_chain_getters_and_setters():
    """Testing the chain accession methods."""
    test_res = Residue(1, "A", list())
    assert test_res.chain is None
    assert test_res.is_missing_chain()

    test_res.chain = "A"
    assert test_res.chain == "A"
    assert not test_res.is_missing_chain()


def test_num_atoms():
    """Making sure the Residue.num_atoms method works correctly."""
    empty_residue = Residue(1, "A", list())
    assert not empty_residue.num_atoms
    assert RESIDUES[0].num_atoms == 16
    res_cpy = deepcopy(RESIDUES[0])
    res_cpy._atoms = []
    assert not empty_residue.num_atoms


def test_renumber_atoms():
    """Ensuring the Residue.renumber_atoms() method works and returns the correct offset number."""

    res_cpy = deepcopy(RESIDUES[0])
    test_idx_list = list(range(1,21))
    res_cpy.renumber_atoms(test_idx_list)
    assert res_cpy.atom_idx_list == list(range(1,17))

    res_cpy = deepcopy(RESIDUES[1])
    res_cpy.renumber_atoms(test_idx_list)
    assert res_cpy.atom_idx_list == list(range(1,15))

def test_remove_atoms_not_in_list():
    """test function works as expected"""
    test_residue = copy.deepcopy(RESIDUES[1])
    test_keep_list = ["C", "H", "CA", "N", "O"]
    assert test_residue.num_atoms == 14
    test_residue.remove_atoms_not_in_list(test_keep_list)
    assert test_residue.num_atoms == 5
    assert set(test_residue.atom_name_list) == set(test_keep_list)

def test_is_connected():
    """test function works as expected"""
    test_residue = copy.deepcopy(RESIDUES[1])
    test_residue.parent = RESIDUES[2].parent
    init_connectivity(test_residue)
    assert test_residue.is_connected()
