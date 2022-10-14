"""Testing the enzy_htp.structure.Residue() class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import itertools
import os
import pytest
import pandas as pd
from typing import List
from copy import deepcopy
from collections import defaultdict
from biopandas.pdb import PandasPdb
from enzy_htp.core.exception import ResidueDontHaveAtom

from enzy_htp.structure import Atom, Residue
from enzy_htp.structure.structure_io import PDBParser
import enzy_htp.structure.structure_operation as so

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"

TEST_PDB_FILE = f"{CURR_DIR}/data/3NIR.pdb"
RESIDUES = so.remove_solvent(PDBParser().get_structure(TEST_PDB_FILE)).residues


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
        atom: Atom = Atom(**row)
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
    r1 = Residue(1, "A", list())
    r2 = Residue(1, "A", list())
    assert r1.is_sequence_eq(r2)
    assert r1.is_sequence_eq(deepcopy(r1))
    r3 = Residue(1, "A", list())
    r3.chain = "B"
    assert r1.is_sequence_eq(r3)
    r3.chain = "A"
    assert r1.is_sequence_eq(r3)
    r3.name = "B"
    assert not r1.is_sequence_eq(r3)


def test_sequence_equivalent_return_false():
    """Testing that the Residue.is_sequence_equivalent() method returns False when desired."""
    r1 = Residue(1, "A", list())
    r2 = Residue(3, "B", list())
    assert not r1.is_sequence_eq(r2)
    r3 = Residue(1, "C", list())
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
    assert res_cpy.renumber_atoms(1) == 16
    assert [aa.atom_number for aa in res_cpy._atoms] == list(range(1, 17))

    res_cpy = deepcopy(RESIDUES[1])
    assert res_cpy.renumber_atoms(10) == 23
    assert [aa.atom_number for aa in res_cpy._atoms] == list(range(10, 24))


def test_renumber_atoms_bad_input():
    """Ensuring the Residue.renumber_atoms() method fails when an invalid start index (<= 0) is given."""

    with pytest.raises(SystemExit) as exe:
        RESIDUES[0].renumber_atoms(0)

    assert exe.type == SystemExit
    assert exe.value.code == 1

    with pytest.raises(SystemExit) as exe:
        RESIDUES[0].renumber_atoms(-1)

    assert exe.type == SystemExit
    assert exe.value.code == 1


def test_deepcopy():
    """Making sure the Residue.__deepcopy__() method returns a deepcopy of the current Residue."""
    res: Residue = deepcopy(RESIDUES[0])
    res_cpy: Residue = deepcopy(res)

    assert id(res) != id(res_cpy)
    for a1, a2 in zip(res.atoms, res_cpy.atoms):
        assert id(a1) != id(a2)
