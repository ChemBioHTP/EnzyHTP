"""Testing the enzy_htp.structure.Atom class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-19
"""
import itertools
import os
import numpy as np
import pandas as pd
from typing import List
from copy import deepcopy
import pytest
from biopandas.pdb import PandasPdb

from enzy_htp.core.file_system import lines_from_file
from enzy_htp.structure.atom import Atom
from enzy_htp.structure.structure_io import PDBParser

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"

def test_deepcopy():
    """test the hehavior of copy.deepcopy on Atom() under a Structure()
    context"""
    stru = PDBParser().get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    atom_list = stru[0][0][1:3] + stru[0][1][1:2] # target for deepcopy
    new_list = deepcopy(atom_list)
    # ensure the list is new
    assert id(atom_list) != id(new_list)
    # ensure every atom is new
    for i, j in zip(atom_list, new_list):
        assert id(i) != id(j)
    # ensure parent of every atom is None
    assert all(i.parent is None for i in new_list)

def test_no_optional_data():
    """test the case pdb dont have atom number"""
    stru = PDBParser().get_structure(f"{DATA_DIR}two_chain_no_optional_data.pdb")
    atom = stru["A"][0][0]
    assert atom.element == "N"
    assert atom.idx is None
    assert atom.b_factor is None
    assert atom.charge is None


def test_element_canonical():
    """test get atom element for C from 1NVG"""
    stru = PDBParser().get_structure(f"{DATA_DIR}1NVG.pdb")
    atom = stru["A"][0][0]
    assert atom.element == "N"

def test_element_non_canonical():
    """test get atom element for first atom of the ligand 4CO from 1Q4T"""
    stru = PDBParser().get_structure(f"{DATA_DIR}1Q4T_ligand_test.pdb")
    atom = stru.ligands[0][0]
    assert atom.element == "N"

def test_element_metal():
    """test get atom element for Zn from 1NVG"""
    stru = PDBParser().get_structure(f"{DATA_DIR}1NVG.pdb")
    atom = stru.metalcenters[0][0]
    assert atom.element == "Zn"

def test_radius():
    stru = PDBParser().get_structure(f"{DATA_DIR}1NVG.pdb")
    atom1 = stru.metalcenters[0][0] #Zn
    atom2 = stru["A"][0][0] #N
    assert atom1.radius() == 0.88
    assert atom2.radius() == 1.32

def test_get_connect_no_h(caplog):
    """get connect with missing H atoms"""
    stru = PDBParser().get_structure(f"{DATA_DIR}1NVG.pdb")
    atom1 = stru["A"][0][0] #N (Nter)
    atom2 = stru["A"][1][0] #N
    atom3 = stru["A"][-1].find_atom_name("C") #C (Cter)
    assert len(atom1.get_connect()) == 1
    assert len(atom2.get_connect()) == 2
    assert len(atom3.get_connect()) == 3
    assert "missing connecting atom " in caplog.text

def test_get_connect_w_h():
    """get connect with no missing H atoms"""
    stru = PDBParser().get_structure(f"{DATA_DIR}1Q4T_peptide_protonated.pdb")
    atom1 = stru["A"][0][0] #N (Nter)
    atom2 = stru["A"][1][0] #N
    atom3 = stru["A"][-1].find_atom_name("C")#N (Cter)
    assert len(atom1.get_connect()) == 4
    assert len(atom2.get_connect()) == 3
    assert len(atom3.get_connect()) == 3

def test_attached_protons():
    """test if the function gives correct attached protons"""
    stru = PDBParser().get_structure(f"{DATA_DIR}1Q4T_peptide_protonated.pdb")
    atom1 = stru["A"][0][0] #N (Nter)
    assert list(map(lambda a: a.name, atom1.attached_protons())) == ['H1', 'H2', 'H3']

# TODO recover when need
@pytest.mark.TODO
def test_distance():
    """Making sure that the .distance() method works."""
    TEST_PDB_NAME = f"{CURR_DIR}/data/atom_test_pdb.pdb"
    RAW_ATOMS: pd.DataFrame = PandasPdb().read_pdb(TEST_PDB_NAME).df["ATOM"]
    atoms: List[Atom] = list(map(lambda pr: Atom(**pr[1]), RAW_ATOMS.iterrows()))
    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=0, z_coord=1)
    assert np.isclose(atom1.distance(atom2), 1)
    assert np.isclose(atom2.distance(atom1), 1)

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=1, z_coord=0)
    assert np.isclose(atom1.distance(atom2), 1)
    assert np.isclose(atom2.distance(atom1), 1)

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=1, y_coord=0, z_coord=0)
    assert np.isclose(atom1.distance(atom2), 1)
    assert np.isclose(atom2.distance(atom1), 1)

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=-1, y_coord=0, z_coord=0)
    assert np.isclose(atom1.distance(atom2), 1)
    assert np.isclose(atom2.distance(atom1), 1)

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=-1, z_coord=0)
    assert np.isclose(atom1.distance(atom2), 1)
    assert np.isclose(atom2.distance(atom1), 1)

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=0, z_coord=-1)
    assert np.isclose(atom1.distance(atom2), 1)
    assert np.isclose(atom2.distance(atom1), 1)

    atom1 = Atom(x_coord=0, y_coord=0, z_coord=0)
    atom2 = Atom(x_coord=0, y_coord=0, z_coord=0)
    assert np.isclose(atom1.distance(atom2), 0)
    assert np.isclose(atom2.distance(atom1), 0)
