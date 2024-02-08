"""Testing the enzy_htp.structure.structure_constraint.xml_io.py submodule
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-01-22
"""
import os

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
import enzy_htp.structure.structure_constraint as stru_cons

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

import pytest


tri_alanine = sp.get_structure(f"{DATA_DIR}/tri_alanine.pdb")

def test_serialized_csts():
    """Function that checks a variety of functionality for a .xml file that has multiple constraints."""
   
    csts = stru_cons.structure_constraints_from_xml(tri_alanine, f"{DATA_DIR}/cst_list.xml")

    assert len(csts) == 6

    assert csts[0].is_residue_pair_constraint()

    assert csts[1].is_distance_constraint()

    assert csts[2].is_angle_constraint()

    assert csts[3].is_dihedral_constraint()

    assert csts[4].is_cartesian_freeze()

    assert csts[5].is_backbone_freeze()
