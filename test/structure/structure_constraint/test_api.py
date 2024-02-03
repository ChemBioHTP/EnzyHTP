"""Testing the enzy_htp.structure.structure_constraint.api.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-25
"""
import pytest
import os
import numpy as np

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
import enzy_htp.structure.structure_constraint as stru_cons

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

def test_angle_current_geom():
    "answer verified using PyMol"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons = stru_cons.create_angle_constraint(
        "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru)
    assert np.isclose(test_cons.current_geometry(), 138.39954225812696, atol=0.0001)

def test_dihedral_current_geom():
    "answer verified using PyMol"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons = stru_cons.create_dihedral_constraint(
        "B.254.CAE", "B.254.H2", "A.101.OE2", "A.101.CA", 0.0, test_stru)
    assert test_cons.current_geometry() == 4.740006673137136

def test_distance_current_geom():
    "answer verified using PyMol"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_cons = stru_cons.create_distance_constraint(
        "B.254.H2", "A.101.OE2", 2.4, test_stru)
    assert test_cons.current_geometry() == 2.0239901185529554
