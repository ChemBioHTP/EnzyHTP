"""Testing the enzy_htp.structure.Atom class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-19
"""
import os
import numpy as np
import pandas as pd
from typing import List
from copy import deepcopy
from biopandas.pdb import PandasPdb

from enzy_htp.core.file_system import lines_from_file
from enzy_htp.structure.atom import Atom

CURR_FILE  = os.path.dirname( os.path.abspath( __file__ ))
TEST_PDB_NAME =f"{CURR_FILE}/data/atom_test_pdb.pdb" #@shaoqz: for the pdb set up maybe we can import from a common test script?
RAW_ATOMS : pd.DataFrame = PandasPdb().read_pdb(TEST_PDB_NAME  ).df['ATOM']
atoms : List[Atom] = list(map(lambda pr: Atom(**pr[1]), RAW_ATOMS.iterrows())) #@shaoqz: what does the **pr mean?



def test_residue_key():
    """Checking that the .residue_key() method works."""
    atom : Atom  = deepcopy(atoms[0])
    assert atom.residue_key() == ".ASP.1"
    atom['chain_id'] = "FIRST"
    atom['residue_name'] = "SECOND"
    atom['residue_number'] = "THIRD"
    assert atom.residue_key() == "FIRST.SECOND.THIRD"

def test_distance():
    """Making sure that the .distance() method works."""
    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=0, z_coord=1)
    assert np.isclose( atom1.distance(atom2), 1 ) 
    assert np.isclose( atom2.distance(atom1), 1 ) 

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=1, z_coord=0)
    assert np.isclose( atom1.distance(atom2), 1 ) 
    assert np.isclose( atom2.distance(atom1), 1 ) 

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=1, y_coord=0, z_coord=0)
    assert np.isclose( atom1.distance(atom2), 1 ) 
    assert np.isclose( atom2.distance(atom1), 1 ) 

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=-1, y_coord=0, z_coord=0)
    assert np.isclose( atom1.distance(atom2), 1 ) 
    assert np.isclose( atom2.distance(atom1), 1 ) 

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=-1, z_coord=0)
    assert np.isclose( atom1.distance(atom2), 1 ) 
    assert np.isclose( atom2.distance(atom1), 1 ) 

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=0, z_coord=-1)
    assert np.isclose( atom1.distance(atom2), 1 ) 
    assert np.isclose( atom2.distance(atom1), 1 ) 

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=0, z_coord=0)
    assert np.isclose( atom1.distance(atom2), 0 ) 
    assert np.isclose( atom2.distance(atom1), 0 ) 

def test_getters_and_setters():
    """Making sure the Atom() object's index operations fully work."""

    atom : Atom = Atom()

    assert not atom['x_coord']
    atom['x_coord'] = 5

    assert np.isclose(atom['x_coord'],5)

    assert not atom['dne'] 

def test_to_pdb_line():
    """Ensuring the .to_pdb_line() method functions properly"""
    actual_lines = list(map(lambda aa: aa.to_pdb_line(), atoms))
    target_lines = lines_from_file( TEST_PDB_NAME )
    for a, l in zip(actual_lines, target_lines):
        print(f'"{a}"')
        print(f'"{l}"')
        assert a == l
    assert actual_lines == target_lines

