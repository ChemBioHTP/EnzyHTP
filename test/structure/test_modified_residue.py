"""Testing the ModifiedResidue() class in enzy_htp.structure. 
Author: Sebastian Stull (sebastian.l.stull@vanderbilt.edu)
Date: 2025-02-18
"""

import os
from enzy_htp.structure import (
    ModifiedResidue,
    Atom,
    )
from enzy_htp.structure.structure_io.pdb_io import PDBParser

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"

def test_find_mainchain():

    a1 = Atom.from_biopandas({"x_coord": 0, "y_coord": 0, "z_coord": 0, "atom_name": "N"})
    a2 = Atom.from_biopandas({"x_coord": 0.96, "y_coord": 0, "z_coord": 0, "atom_name": "CA"})
    a3 = Atom.from_biopandas({"x_coord": -0.24, "y_coord": 0.93, "z_coord": 1, "atom_name": "C"})
    a4 = Atom.from_biopandas({"x_coord": -0.26, "y_coord": 1.94, "z_coord": 1, "atom_name": "CB"})
    a5 = Atom.from_biopandas({"x_coord": 0.26, "y_coord": -1.94, "z_coord": -1, "atom_name": "O"})

    # create fake connectivity
    a1.connect = [(a2, "s")]
    a2.connect = [(a1, "s"), (a3, "s")]
    a3.connect = [(a4, "s")]
    a4.connect = [(a5, "s")]
    a5.connect = [(a4, "s")]

    ncaa_1 = ModifiedResidue(
        residue_idx=1,
        residue_name="TST",
        atoms=[a1, a2, a3, a4, a5],
    )

    # assert main chain is correct
    mc = ncaa_1.find_mainchain()
    res = [aa.name for aa in mc]

    assert res == ["N", "CA", "C"]


def test_modified_residue_clone_with_connectivity():
    """Test ModifiedResidue.clone() preserves connectivity and all attributes"""
    from enzy_htp.structure.atom import Atom
    
    # Create atoms with connectivity
    atom1 = Atom("N", (0.0, 0.0, 0.0))
    atom2 = Atom("CA", (1.0, 0.0, 0.0))
    atom3 = Atom("C", (2.0, 0.0, 0.0))
    atom4 = Atom("CB", (1.0, 1.0, 0.0))
    
    # Set up connectivity
    atom1.connect_to(atom2, "single")
    atom2.connect_to(atom3, "single")
    atom2.connect_to(atom4, "single")
    
    # Create ModifiedResidue with attributes
    mod_res = ModifiedResidue(1, "MSE", [atom1, atom2, atom3, atom4], 
                             net_charge=0, multiplicity=1, mainchain_atom_names=["N", "CA", "C"])
    
    # Clone with connectivity
    mod_res_clone = mod_res.clone(with_connectivity=True)
    
    # Verify basic properties
    assert isinstance(mod_res_clone, ModifiedResidue)
    assert mod_res_clone.idx == mod_res.idx
    assert mod_res_clone.name == mod_res.name
    assert len(mod_res_clone.atoms) == len(mod_res.atoms)
    
    # Verify attributes are preserved
    assert mod_res_clone.net_charge == mod_res.net_charge
    assert mod_res_clone.multiplicity == mod_res.multiplicity
    
    # Verify mainchain atoms are preserved
    assert len(mod_res_clone.mainchain_atoms) == 3
    mainchain_names = [atom.name for atom in mod_res_clone.mainchain_atoms]
    assert mainchain_names == ["N", "CA", "C"]
    
    # Verify connectivity is preserved
    cloned_atoms = {atom.name: atom for atom in mod_res_clone.atoms}
    assert cloned_atoms["N"].is_connected()
    assert len(cloned_atoms["N"].connect) == 1
    assert cloned_atoms["CA"].is_connected()
    assert len(cloned_atoms["CA"].connect) == 3  # Connected to N, C, and CB
    
    # Check specific connections
    ca_connections = {atom.name: bond for atom, bond in cloned_atoms["CA"].connect}
    assert "N" in ca_connections and ca_connections["N"] == "single"
    assert "C" in ca_connections and ca_connections["C"] == "single"
    assert "CB" in ca_connections and ca_connections["CB"] == "single"


def test_modified_residue_clone_without_connectivity():
    """Test ModifiedResidue.clone() without connectivity preservation"""
    from enzy_htp.structure.atom import Atom
    
    # Create atoms with connectivity
    atom1 = Atom("N", (0.0, 0.0, 0.0))
    atom2 = Atom("CA", (1.0, 0.0, 0.0))
    atom1.connect_to(atom2, "single")
    
    # Create ModifiedResidue
    mod_res = ModifiedResidue(1, "MSE", [atom1, atom2], net_charge=-1, multiplicity=2)
    
    # Clone without connectivity
    mod_res_clone = mod_res.clone(with_connectivity=False)
    
    # Verify basic properties
    assert isinstance(mod_res_clone, ModifiedResidue)
    assert mod_res_clone.net_charge == mod_res.net_charge
    assert mod_res_clone.multiplicity == mod_res.multiplicity
    
    # Verify connectivity is not preserved
    for atom in mod_res_clone.atoms:
        assert not atom.is_connected()