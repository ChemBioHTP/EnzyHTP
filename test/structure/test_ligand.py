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

# TODO recover these test
@pytest.mark.TODO
def test_net_charge():
    """Checking that the net_charge attribute works properly for both default and set values."""
    lig = Ligand("A.B.10", list())
    assert lig.get_net_charge() is None
    lig = Ligand("A.B.10", list(), net_charge=1.0)
    assert np.isclose(lig.get_net_charge(), 1.0)

@pytest.mark.TODO
def test_build_bad_fname():
    """Ensuring the Ligand.build() method fails when a non .pdb files is supplied."""
    lig = Ligand("A.B.10", list())
    with pytest.raises(SystemExit) as exe:
        lig.build("test.txt")

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1

@pytest.mark.TODO
def test_clone():
    """Checking that the Ligand.clone() method returns a deepcopy of the current Ligand()."""
    dummy_atoms = [Atom(line_idx=1), Atom(line_idx=2)]
    lig = Ligand("A.B.10", dummy_atoms, net_charge=1.0)
    lig_cpy = lig.clone()
    assert isinstance(lig_cpy, Ligand)
    assert id(lig) != id(lig_cpy)
    for a1, a2 in zip(lig.atoms, lig_cpy.atoms):
        assert id(a1) != id(a2)
