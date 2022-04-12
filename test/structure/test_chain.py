"""Testing the functionality of the enzy_htp.structure.Chain class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-20
"""
import os

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"


from enzy_htp.structure import Chain, Residue, structure_from_pdb


def test_proper_ctor_behavior():
    """Making sure that the default Chain() works."""
    chain = Chain("test", [])
    assert chain.name() == "test"
    assert chain.empty()
    assert not chain.residues()


# TODO(CJ): need to add tests for checking if the chain has ceratin residues/residue types.


def test_same_sequence_equal():
    """Ensuring the Chain.same_sequence() method works for equivalent chains."""
    residues = [
        Residue("A.A.1", list()),
        Residue("A.A.2", list()),
        Residue("A.B.3", list()),
        Residue("A.D.5", list()),
    ]
    chain1 = Chain("A", residues)
    chain2 = Chain("A", residues)
    assert chain1.same_sequence(chain2)
    assert chain2.same_sequence(chain1)

    residues_alt = [
        Residue("B.A.1", list()),
        Residue("B.A.2", list()),
        Residue("B.B.3", list()),
        Residue("B.D.5", list()),
    ]
    chain3 = Chain("A", residues_alt)
    assert chain3.same_sequence(chain2)
    assert chain2.same_sequence(chain3)


def test_same_sequence_not_equal():
    """Ensuring the Chain.same_sequence() method works for non-equivalent chains."""
    residues = [
        Residue("A.A.1", list()),
        Residue("A.A.2", list()),
        Residue("A.B.3", list()),
        Residue("A.D.5", list()),
    ]

    chain1 = Chain("A", residues)
    chain2 = Chain("B", residues)
    assert not chain1.same_sequence(chain2)
    assert not chain2.same_sequence(chain1)

    residues_alt1 = [
        Residue("A.B.1", list()),
        Residue("A.A.2", list()),
        Residue("A.B.3", list()),
        Residue("A.D.5", list()),
    ]
    chain3 = Chain("A", residues_alt1)
    assert not chain3.same_sequence(chain1)
    assert not chain3.same_sequence(chain1)

    residues_alt2 = [
        Residue("A.A.3", list()),
        Residue("A.A.2", list()),
        Residue("A.B.3", list()),
        Residue("A.D.5", list()),
    ]
    chain4 = Chain("A", residues_alt2)
    assert not chain4.same_sequence(chain1)
    assert not chain1.same_sequence(chain4)

    residues.pop()
    residues.pop()
    chain5 = Chain("A", residues)
    assert not chain5.same_sequence(chain1)
    assert not chain1.same_sequence(chain5)


def test_rename():
    """Testing the Chain.rename() method that gives the Chain() a new name."""
    residues = [
        Residue("A.A.1", list()),
        Residue("A.A.2", list()),
        Residue("A.B.3", list()),
        Residue("A.D.5", list()),
    ]
    chain = Chain("A", residues)

    for res in chain.residues():
        assert res.chain() == "A"

    assert chain.name() == "A"
    chain.rename("B")

    assert chain.name() == "B"

    for res in chain.residues():
        assert res.chain() == "B"


def test___getitem__():
    """Testing integer __getitem__() indexing for the Chain() object."""
    res_keys = ["A.A.1", "A.A.2", "A.B.3", "A.D.5"]
    residues = list(map(lambda k: Residue(k, list()), res_keys))
    chain = Chain("A", residues)

    for ridx, rkey in enumerate(res_keys):
        assert chain[ridx].residue_key == rkey


def test___delitem__():
    """Testing integer __delitem__() indexing for the Chain() object."""
    res_keys = ["A.A.1", "A.A.2", "A.B.3", "A.D.5"]
    residues = list(map(lambda k: Residue(k, list()), res_keys))
    chain = Chain("A", residues)

    del chain[2]
    del res_keys[2]

    for ridx, rkey in enumerate(res_keys):
        assert chain[ridx].residue_key == rkey


def test_num_atoms():
    """Finds the total number of Atom() objects in a given Chain() object."""
    res_keys = ["A.A.1", "A.A.2", "A.B.3", "A.D.5"]
    residues = list(map(lambda k: Residue(k, list()), res_keys))
    empty_chain = Chain("A", residues)
    assert empty_chain.num_atoms() == 0


def test_is_metal():
    """Checks if the Chain.is_metal() returns correct answers for both True and False cases."""
    pdb_file = f"{DATA_DIR}/1NVG.pdb"
    structure: Structure = structure_from_pdb(pdb_file)
    chain: Chain = structure.chains()[0]
    assert chain.is_metal()
    del chain[-1]
    del chain[-1]
    assert not chain.is_metal()
