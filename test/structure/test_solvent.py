"""Testing the Solvent() class that represents a Solvent() in a .pdb file. Derived from the Residue() class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-04-05
"""
import pytest
from enzy_htp.chemical import enum as renum
from enzy_htp.structure import Solvent, residue_to_solvent, Atom, Residue


def test_constant_data():
    """Ensuring the constant data for the Solvent() class is encoded properly."""
    solvent = Solvent(1, "HOH", list())
    assert solvent.is_solvent()
    assert not solvent.is_metal()
    assert not solvent.is_ligand()
    assert not solvent.is_canonical()
    assert solvent.rtype == renum.ResidueType.SOLVENT


def test_residue_to_solvent():
    """Making sure the enzy_htp.structure.solvent.residue_to_solvent() method works."""
    residue = Residue(1, "HOH", list())
    solvent: Solvent = residue_to_solvent(residue)
    assert isinstance(solvent, Solvent)
    assert id(residue) != id(solvent)
    for a1, a2 in zip(solvent.atoms, residue.atoms):
        assert id(a1) != id(a2)
