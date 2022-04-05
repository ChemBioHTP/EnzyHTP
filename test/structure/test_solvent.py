"""Testing the Solvent() class that represents a Solvent() in a .pdb file. Derived from the Residue() class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-04-05
"""
from enzy_htp.chemical import enum as renum
from enzy_htp.structure import Solvent, residue_to_solvent, Atom, Residue


def test_constant_data():
    """Ensuring the constant data for the Solvent() class is encoded properly."""
    solvent = Solvent("A.HOH.1", list())
    assert solvent.is_rd_solvent()
    assert not solvent.is_metal()
    assert not solvent.is_ligand()
    assert not solvent.is_canonical()
    assert solvent.rtype() == renum.ResidueType.SOLVENT


def test_clone():
    """Making sure the Solvent.clone() method makes a proper deepcopy."""
    solvent = Solvent("A.HOH.1", [Atom(line_idx=1), Atom(line_idx=2)])
    solvent_cpy: Solvent = solvent.clone()
    assert id(solvent) != id(solvent_cpy)
    for a1, a2 in zip(solvent.atom_list(), solvent_cpy.atom_list()):
        assert id(a1) != id(a2)


def test_residue_to_solvent():
    """Making sure the enzy_htp.structure.solvent.residue_to_solvent() method works."""
    residue = Residue("A.HOH.1", [Atom(line_idx=1), Atom(line_idx=2)])
    solvent: Solvent = residue_to_solvent(residue)
    assert isinstance(solvent, Solvent)
    assert id(residue) != id(solvent)
    for a1, a2 in zip(solvent.atom_list(), residue.atom_list()):
        assert id(a1) != id(a2)
