"""Testing the Mol2Parser class in the enzy_htp.structure.structure_io.mol2_io.py

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-10-23
"""
from pathlib import Path
import tempfile
import re

import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_io import Mol2Parser
from enzy_htp.structure.mol_desc_data import MolDescData
from enzy_htp.structure.structure_io import PrepinParser

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/../data"

def test_mol2_parser_get_stru():
    """make sure function works as expected"""
    test_prepin = f"{DATA_DIR}/TYQ.mol2"
    test_stru = Mol2Parser().get_structure(test_prepin)
    assert len(test_stru.residues) == 1

def test_write_from_mol_desc_data_simple():
    """Test writing a simple MolDescData to a mol2 file."""
    data = MolDescData(
        name="LIG",
        atoms=[
            {'id': 1, 'atom_name': 'C1', 'atom_type': 'C.3', 'charge': -0.1, 'coords': [0.0, 0.0, 0.0]},
            {'id': 2, 'atom_name': 'O1', 'atom_type': 'O.2', 'charge': -0.5, 'coords': [1.4, 0.0, 0.0]},
        ],
        bonds=[
            {'atom1_id': 1, 'atom2_id': 2, 'bond_type': '1'},
        ]
    )
    
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.mol2') as tmp:
        Mol2Parser.write_from_mol_desc_data(data, tmp.name)
        tmp.seek(0)
        content = tmp.read()
        print(content)  # For debugging purposes
        assert "@<TRIPOS>MOLECULE" in content
        assert "LIG" in content
        assert "    2     1     1     0     0" in content
        assert "@<TRIPOS>ATOM" in content
        assert "      1 C1           0.0000     0.0000     0.0000 C.3        1 LIG         -0.1000" in content
        assert "      2 O1           1.4000     0.0000     0.0000 O.2        1 LIG         -0.5000" in content
        assert "@<TRIPOS>BOND" in content
        assert re.search(r"1\s+1\s+2\s+1", content)
        assert "@<TRIPOS>SUBSTRUCTURE" in content

def test_write_from_mol_desc_data_and_read_back():
    """Test writing and reading back a molecule."""
    test_prepin = f"{DATA_DIR}/ligand_H5J.prepin"
    mol_desc_from_prepin = PrepinParser.get_mol_desc_data(test_prepin)
    
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.mol2') as tmp:
        Mol2Parser.write_from_mol_desc_data(mol_desc_from_prepin, tmp.name)
        
        # Now read the mol2 file back and check consistency
        mol_desc_from_mol2 = Mol2Parser.get_ligand(tmp.name)
        
        assert mol_desc_from_prepin.name == mol_desc_from_mol2.name[:3] # mol2 parser takes first 3 letters
        assert len(mol_desc_from_prepin.atoms) == len(mol_desc_from_mol2.atoms)
        
        # Check atom names and charges (mol2 format might have different precision)
        for atom1 in mol_desc_from_prepin.atoms:
            atom2 = next((a for a in mol_desc_from_mol2.atoms if a.name == atom1['atom_name']), None)
            assert atom2 is not None
            assert round(atom1['charge'], 3) == round(atom2.charge, 3)