"""Testing the enzy_htp.structure.Atom class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-19
"""
import os
import numpy as np
import pandas as pd
from typing import List
from copy import deepcopy
import pytest
from biopandas.pdb import PandasPdb

from enzy_htp.core.file_system import lines_from_file
from enzy_htp.structure.atom import Atom

CURR_FILE = os.path.dirname(os.path.abspath(__file__))
TEST_PDB_NAME = f"{CURR_FILE}/data/atom_test_pdb.pdb"
RAW_ATOMS: pd.DataFrame = PandasPdb().read_pdb(TEST_PDB_NAME).df["ATOM"]
atoms: List[Atom] = list(map(lambda pr: Atom(**pr[1]), RAW_ATOMS.iterrows()))


# TODO recover when need
@pytest.mark.TODO
def test_distance():
    """Making sure that the .distance() method works."""
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
