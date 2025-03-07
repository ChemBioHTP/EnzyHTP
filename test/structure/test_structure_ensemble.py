"""Testing the StructureEnsemble() class. 
Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-01-03
"""
import os
import pytest
import numpy as np
from copy import deepcopy
import enzy_htp.chemical as chem
from enzy_htp.core import file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate, get_itself
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp.structure import (
    Structure,
    Residue,
    Chain,
    Atom,
    Ligand,
)
from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.structure.structure_io import PrmtopParser
from enzy_htp._interface.amber_interface import AmberMDCRDParser, AmberNCParser

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
WORK_DIR = f"{CURR_DIR}/work_dir/"
DATA_DIR = f"{CURR_DIR}/data/"
sp = PDBParser()


def test_structures():
    """test the function in name
    use Structure and Amber as example
    check using manually confirmed answer"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_stru.pbc_box_shape = (66.957, 66.957, 66.957, 109.471219, 109.471219, 109.471219)
    test_traj = f"{DATA_DIR}KE_07_R7_2_S_10f.mdcrd"
    test_prmtop = f"{DATA_DIR}KE_07_R7_2_S_10f.prmtop"
    test_esm = StructureEnsemble(
        topology=test_stru,
        top_parser=get_itself,
        coordinate_list=test_traj,
        coord_parser=AmberMDCRDParser(test_prmtop).get_coordinates,
    )
    answer = [
        ([32.586, 55.789, 30.602], [39.661, 27.409, 36.51], (66.957, 66.957, 66.957)),
        ([34.88, 53.201, 22.495], [36.392, 27.392, 38.547], (66.937, 66.937, 66.937)),
        ([37.41, 54.562, 25.501], [35.305, 27.228, 38.266], (66.897, 66.897, 66.897)),
        ([27.814, 51.521, 21.567], [37.699, 30.371, 39.822], (66.918, 66.918, 66.918)),
        ([30.507, 43.798, 13.187], [35.46, 35.527, 40.971], (66.848, 66.848, 66.848)),
        ([30.84, 19.298, 14.47], [35.6, 39.992, 38.233], (66.882, 66.882, 66.882)),
        ([29.9, 19.785, 17.46], [36.833, 40.408, 37.839], (66.855, 66.855, 66.855)),
        ([29.248, 19.402, 17.4], [38.487, 38.876, 37.831], (66.923, 66.923, 66.923)),
        ([25.089, 14.023, 21.701], [36.193, 41.793, 35.419], (66.842, 66.842, 66.842)),
        ([22.332, 17.39, 19.672], [38.004, 39.988, 36.79], (66.896, 66.896, 66.896)),
        ([39.391, 19.22, 17.043], [37.366, 39.457, 37.706], (66.934, 66.934, 66.934)),
    ]

    for stru, (answer_1, answer_m1, answer_pbc_edges) in zip(test_esm.structures(), answer):
        assert stru.atoms[0].coord == answer_1
        assert stru.atoms[-1].coord == answer_m1
        assert stru.pbc_box_shape[:3] == answer_pbc_edges

def test_structures_rm_solvent():
    """test the function in name
    use Structure and Amber as example
    remove_solvent=True
    check using manually confirmed answer"""
    test_prmtop = f"{DATA_DIR}stru_w_solvent.prmtop"
    test_traj = f"{DATA_DIR}traj_w_solvent_3f.nc"
    test_esm = StructureEnsemble(
        topology=test_prmtop,
        top_parser=PrmtopParser({
                "FAD" : (0,1),
                "ACP" : (0,1),
                }
            ).get_structure,
        coordinate_list=test_traj,
        coord_parser=AmberNCParser(test_prmtop).get_coordinates,
    )
    answer = [
        ([69.168, 45.518, 83.438], [33.796, 76.501,  9.040]),
        ([68.969, 46.198, 82.766], [32.082, 77.429, 12.821]),
        ([68.605, 45.599, 82.123], [30.553, 75.828,  9.951]),
    ]

    for stru, (answer_1, answer_m1) in zip(test_esm.structures(), answer):
        assert stru.atoms[0].coord == answer_1
        assert stru.atoms[-1].coord == answer_m1
