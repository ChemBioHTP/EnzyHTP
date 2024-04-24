"""Testing enzy_htp.analysis.spi_metric.py

Author: Chris Jurich <chris.jurich@vanderbilt.edu> 
Date: 2024-02-24
"""

import glob
import pytest
import os
import numpy as np

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_region import create_region_from_selection_pattern
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.analysis import bond_dipole
from enzy_htp import interface
from enzy_htp import PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()


def test_spi_consistent_with_old_enzyhtp():
    print(DATA_DIR)
    assert False
    pass
