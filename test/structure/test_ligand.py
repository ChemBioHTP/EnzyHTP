"""Testing the Ligand() class in enzy_htp.structure. 
Author: Chris Jurich <cjurich2@huskers.unl.edu>
Date: 2022-04-03
"""
import os
import pytest
import numpy as np

from enzy_htp.chemical import enum as renum
from enzy_htp.structure import Ligand, Atom, PDBParser

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"


def test_constat_data():
    """Testing a variety of constant data methods that should work."""
    lig = Ligand(10, "X", list())
    assert lig.is_ligand()
    assert not lig.is_canonical()
    assert not lig.is_metal()
    assert lig.rtype == renum.ResidueType.LIGAND


def test_fix_atom_names():
    """test using a ligand with element symbol as atom names"""
    stru = PDBParser().get_structure(f"{DATA_DIR}just_ligand_badname.pdb")
    ligand = stru.ligands[0]
    ligand.fix_atom_names()
    assert ligand.atom_name_list == ['C', 'F', 'O', 'C1', 'O1', 'H', 'H1']


def test_net_charge():
    """Checking that the net_charge attribute works properly for both default and set values."""
    lig = Ligand(10, "B", list())
    assert lig.net_charge is None
    lig = Ligand(10, "B", list(), net_charge=1.0)
    assert np.isclose(lig.net_charge, 1.0)


def test_clone():
    """Checking that the Ligand.clone() method returns a deepcopy of the current Ligand()."""
    dummy_atoms = [
        Atom({
            'x_coord': 0,
            'y_coord': 0,
            'z_coord': 0,
            'atom_name': 'dummy'
        }),
        Atom({
            'x_coord': 0,
            'y_coord': 0,
            'z_coord': 0,
            'atom_name': 'dummy'
        })
    ]
    lig = Ligand(10, "B", dummy_atoms, net_charge=1.0)
    lig_cpy = lig.clone()
    assert isinstance(lig_cpy, Ligand)
    assert id(lig) != id(lig_cpy)
    for a1, a2 in zip(lig.atoms, lig_cpy.atoms):
        assert id(a1) != id(a2)
