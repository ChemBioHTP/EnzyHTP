"""Testing the enzy_htp.structure.Atom class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: QZ Shao <shaoqz@icloud.com>
Date: 2022-03-19
"""
import itertools
import logging
import os
import enzy_htp
import numpy as np
import pandas as pd
from typing import List
from copy import deepcopy
import pytest
from biopandas.pdb import PandasPdb

from enzy_htp.core.file_system import lines_from_file
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure.atom import Atom
from enzy_htp.structure.structure_io import PDBParser
from enzy_htp.structure.structure_operation import init_connectivity

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"
_LOGGER.propagate = 1

def test_deepcopy():
    """test the hehavior of copy.deepcopy on Atom() under a Structure()
    context"""
    stru = PDBParser().get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    atom_list = stru[0][0][1:3] + stru[0][1][1:2]  # target for deepcopy
    new_list = deepcopy(atom_list)
    # ensure the list is new
    assert id(atom_list) != id(new_list)
    # ensure every atom is new
    for i, j in zip(atom_list, new_list):
        assert id(i) != id(j)
    # ensure parent of every atom is None
    assert all(i.parent is None for i in new_list)

def test_deepcopy_more_than_once():
    """test the hehavior of copy.deepcopy on Atom() under a Structure()
    context and if copied more than once"""
    stru = PDBParser().get_structure(f"{DATA_DIR}12E8_small_four_chain.pdb")
    atom_list = stru[0][0][1:3] + stru[0][1][1:2]  # target for deepcopy
    assert atom_list[0].name != "X"
    new_list = deepcopy(atom_list)
    new_list[0].name = "X"
    new_new_list = deepcopy(new_list)
    assert new_new_list[0].name == "X"
    # ensure parent of every atom is None (so not using the python default deepcopy)
    assert all(i.parent is None for i in new_list)
    assert all(i.parent is None for i in new_new_list)

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
    atom1 = stru.metalcenters[0][0]  #Zn
    atom2 = stru["A"][0][0]  #N
    assert atom1.radius() == 0.88
    assert atom2.radius() == 1.32


def test_attached_protons():
    """test if the function gives correct attached protons"""
    stru = PDBParser().get_structure(f"{DATA_DIR}1Q4T_peptide_protonated.pdb")
    atom1: Atom = stru["A"][0][0]  #N (Nter)
    init_connectivity(atom1)
    assert list(map(lambda a: a.name, atom1.attached_protons())) == ['H1', 'H2', 'H3']


def test_distance_to():
    """Making sure that the .distance_to() method works."""
    atom1 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 0, 'atom_name': 'DUMMY'})

    atom2 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 1, 'atom_name': 'DUMMY'})
    assert np.isclose(atom1.distance_to(atom2), 1)
    assert np.isclose(atom2.distance_to(atom1), 1)

    atom2 = Atom({'x_coord': 0, 'y_coord': 1, 'z_coord': 0, 'atom_name': 'DUMMY'})
    assert np.isclose(atom1.distance_to(atom2), 1)
    assert np.isclose(atom2.distance_to(atom1), 1)

    atom2 = Atom({'x_coord': 1, 'y_coord': 0, 'z_coord': 0, 'atom_name': 'DUMMY'})
    assert np.isclose(atom1.distance_to(atom2), 1)
    assert np.isclose(atom2.distance_to(atom1), 1)

    atom2 = Atom({'x_coord': -1, 'y_coord': 0, 'z_coord': 0, 'atom_name': 'DUMMY'})
    assert np.isclose(atom1.distance_to(atom2), 1)
    assert np.isclose(atom2.distance_to(atom1), 1)

    atom2 = Atom({'x_coord': 0, 'y_coord': -1, 'z_coord': 0, 'atom_name': 'DUMMY'})
    assert np.isclose(atom1.distance_to(atom2), 1)
    assert np.isclose(atom2.distance_to(atom1), 1)

    atom2 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': -1, 'atom_name': 'DUMMY'})
    assert np.isclose(atom1.distance_to(atom2), 1)
    assert np.isclose(atom2.distance_to(atom1), 1)

    atom2 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 0, 'atom_name': 'DUMMY'})
    assert np.isclose(atom1.distance_to(atom2), 0)
    assert np.isclose(atom2.distance_to(atom1), 0)

def test_check_connect_setter_data_type_correct():
    """check using correct data type"""
    atom1 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 0, 'atom_name': 'DUMMY'})
    atom2 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 1, 'atom_name': 'DUMMY'})

    atom1.connect = [(atom2, "s")]
    atom1.connect = [(atom2, None)]

def test_check_connect_setter_data_type_wrong():
    """check using several wrong data type"""
    atom1 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 0, 'atom_name': 'DUMMY'})
    atom2 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 1, 'atom_name': 'DUMMY'})

    with pytest.raises(TypeError) as exe:
        atom1.connect = [[atom2, "s"]]
    assert exe.value.args[0] == "Assigning wrong data type for connect. Correct data type: [[Atom(), 'bond_order_info'], ...]"

    with pytest.raises(TypeError):
        atom1.connect = [atom2, "s"]

    with pytest.raises(TypeError):
        atom1.connect = [atom2]

    with pytest.raises(TypeError):
        atom1.connect = [(atom2, "s"), atom2]

def test_connect_to(caplog):
    """test using example"""
    _LOGGER.setLevel(logging.DEBUG)
    atom1 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 0, 'atom_name': 'DUMMY'})
    atom2 = Atom({'x_coord': 0, 'y_coord': 0, 'z_coord': 1, 'atom_name': 'DUMMY'})

    atom1.connect_to(atom2)
    atom1.connect_to(atom2)
    assert "already in connect" in caplog.text

    assert atom2 in atom1.connect_atoms
    assert atom1 in atom2.connect_atoms
