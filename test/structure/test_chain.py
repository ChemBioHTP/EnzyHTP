"""Testing the functionality of the enzy_htp.structure.Chain class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-20
"""
from copy import deepcopy
import logging
import os
import pytest
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Chain, Residue, PDBParser, Structure

_LOGGER.setLevel(logging.DEBUG)
CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"
sp = PDBParser()


def test_sequence():
    """test getting the sequence of the chain"""
    pdb_file_path = f"{DATA_DIR}two_chain.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)

    assert stru[0].sequence == "AT"
    assert stru[1].sequence == "GGNLP"


def test_sequence_noncanonical():
    """test getting the sequence of the chain"""
    pdb_file_path = f"{DATA_DIR}5JT3_noncanonical_test.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    assert stru[0].sequence == "KRMLNTGYSLNNVHIDYVPTV TPO A"


def test_is_same_coord():
    """Testing that the Chain.is_same_coord() method workds correctly"""
    pdb_file_path = f"{DATA_DIR}two_chain.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    ch1, ch2 = deepcopy(stru.chains[0]), deepcopy(stru.chains[0])
    assert ch1.is_same_coord(ch2)
    assert ch2.is_same_coord(ch1)
    ch2.atoms[0].coord = (ch2.atoms[0].coord[0] + 0.05, ch2.atoms[0].coord[1],
                          ch2.atoms[0].coord[2])
    assert not ch1.is_same_coord(ch2)
    assert ch1.is_same_coord(ch2, 0.10)


def test_proper_ctor_behavior():
    """Making sure that the default Chain() works."""
    chain = Chain("test", [])
    assert chain.name == "test"
    assert chain.is_empty()
    assert not chain.residues


# TODO(CJ): need to add tests for checking if the chain has ceratin residues/residue types.


def test_same_sequence_equal():
    """Ensuring the Chain.is_same_sequence() method works for equivalent chains."""
    residues = [
        Residue(1, "A", list()),
        Residue(2, "A", list()),
        Residue(3, "B", list()),
        Residue(5, "D", list()),
    ]
    chain1 = Chain("A", residues)
    chain2 = Chain("A", residues)
    assert chain1.is_same_sequence(chain2)
    assert chain2.is_same_sequence(chain1)

    residues_alt = [
        Residue(1, "A", list()),
        Residue(2, "A", list()),
        Residue(3, "B", list()),
        Residue(5, "D", list()),
    ]
    chain3 = Chain("A", residues_alt)
    assert chain3.is_same_sequence(chain2)
    assert chain2.is_same_sequence(chain3)


@pytest.mark.TODO
def test_same_sequence_not_equal():
    """Ensuring the Chain.is_same_sequence() method works for non-equivalent chains."""
    residues = [
        Residue(1, "A", list()),
        Residue(2, "A", list()),
        Residue(3, "B", list()),
        Residue(5, "D", list()),
    ]

    chain1 = Chain("A", residues)

    residues_alt1 = [
        Residue(1, "B", list()),
        Residue(2, "A", list()),
        Residue(3, "B", list()),
        Residue(5, "D", list()),
    ]
    chain3 = Chain("A", residues_alt1)
    assert not chain3.is_same_sequence(chain1)
    assert not chain3.is_same_sequence(chain1)

    rr = deepcopy(residues)
    rr.pop()
    rr.pop()
    chain5 = Chain("A", rr)
    assert not chain5.is_same_sequence(chain1)
    assert not chain1.is_same_sequence(chain5)


def test_rename():
    """Testing the Chain.rename() method that gives the Chain() a new name."""
    residues = [
        Residue(1, "A", list()),
        Residue(2, "A", list()),
        Residue(3, "B", list()),
        Residue(5, "D", list()),
    ]
    chain = Chain("A", residues)

    for res in chain.residues:
        assert res.chain.name == "A"

    assert chain.name == "A"
    chain.name = "B"

    assert chain.name == "B"

    for res in chain.residues:
        assert res.chain.name == "B"


def test_num_atoms():
    """Finds the total number of Atom() objects in a given Chain() object."""
    res_keys = ["A.A.1", "A.A.2", "A.B.3", "A.D.5"]
    residues = [
        Residue(1, "A", list()),
        Residue(2, "A", list()),
        Residue(3, "B", list()),
        Residue(5, "D", list()),
    ]
    empty_chain = Chain("A", residues)
    assert empty_chain.num_atoms == 0


def test_has_metal():
    """Checks if the Chain.has_metal() returns correct answers for both True and False cases."""
    pdb_file = f"{DATA_DIR}/1NVG.pdb"
    structure: Structure = PDBParser.get_structure(pdb_file)
    assert not structure.chains[0].has_metal()
    assert structure.chains[1].has_metal()
