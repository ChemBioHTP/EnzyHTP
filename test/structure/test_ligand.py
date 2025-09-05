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
        Atom.from_biopandas({
            'x_coord': 0,
            'y_coord': 0,
            'z_coord': 0,
            'atom_name': 'dummy'
        }),
        Atom.from_biopandas({
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


def test_ligand_clone_with_connectivity():
    """Test Ligand.clone() preserves connectivity and all attributes"""
    from enzy_htp.structure.atom import Atom
    
    # Create atoms with connectivity
    atom1 = Atom("C1", (0.0, 0.0, 0.0))
    atom2 = Atom("C2", (1.0, 0.0, 0.0))
    atom3 = Atom("O1", (0.0, 1.0, 0.0))
    
    # Set up connectivity
    atom1.connect_to(atom2, "single")
    atom1.connect_to(atom3, "double")
    
    # Create Ligand with attributes
    lig = Ligand(1, "LIG", [atom1, atom2, atom3], net_charge=-1, multiplicity=2)
    lig.bonds = [{"type": "single"}, {"type": "double"}]
    lig.placement_method = "docking"
    
    # Clone with connectivity
    lig_clone = lig.clone(with_connectivity=True)
    
    # Verify basic properties
    assert isinstance(lig_clone, Ligand)
    assert lig_clone.idx == lig.idx
    assert lig_clone.name == lig.name
    assert len(lig_clone.atoms) == len(lig.atoms)
    
    # Verify attributes are preserved
    assert lig_clone.net_charge == lig.net_charge
    assert lig_clone.multiplicity == lig.multiplicity
    assert lig_clone.bonds == lig.bonds
    assert lig_clone.placement_method == lig.placement_method
    
    # Verify connectivity is preserved
    cloned_atoms = {atom.name: atom for atom in lig_clone.atoms}
    assert cloned_atoms["C1"].is_connected()
    assert len(cloned_atoms["C1"].connect) == 2
    
    # Check specific connections
    c1_connections = {atom.name: bond for atom, bond in cloned_atoms["C1"].connect}
    assert "C2" in c1_connections and c1_connections["C2"] == "single"
    assert "O1" in c1_connections and c1_connections["O1"] == "double"


def test_ligand_clone_without_connectivity():
    """Test Ligand.clone() without connectivity preservation"""
    from enzy_htp.structure.atom import Atom
    
    # Create atoms with connectivity
    atom1 = Atom("C1", (0.0, 0.0, 0.0))
    atom2 = Atom("C2", (1.0, 0.0, 0.0))
    atom1.connect_to(atom2, "single")
    
    # Create Ligand
    lig = Ligand(1, "LIG", [atom1, atom2], net_charge=0, multiplicity=1)
    
    # Clone without connectivity
    lig_clone = lig.clone(with_connectivity=False)
    
    # Verify basic properties
    assert isinstance(lig_clone, Ligand)
    assert lig_clone.net_charge == lig.net_charge
    assert lig_clone.multiplicity == lig.multiplicity
    
    # Verify connectivity is not preserved
    for atom in lig_clone.atoms:
        assert not atom.is_connected()
