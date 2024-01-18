"""Testing enzy_htp._interface.multiwfn_interface.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-01-16
"""
import re
import pytest
import os
import numpy as np

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.chemical.level_of_theory import QMLevelOfTheory, MMLevelOfTheory
from enzy_htp._config.armer_config import ARMerConfig
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp.structure.structure_region import create_region_from_selection_pattern
from enzy_htp.structure.structure_operation.charge import init_charge
from enzy_htp import PDBParser
from enzy_htp import interface
from enzy_htp._interface.gaussian_interface import GaussianSinglePointEngine

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()
mi = interface.multiwfn


def test_parse_two_center_dp_moments():
    """as name. use an example file from another project."""
    test_dipole_file = f"{DATA_DIR}bond_dipole_LMO_multiwfn_example.wfndip"

    test_result = mi.parse_two_center_dp_moments(test_dipole_file)

    assert sum([len(v) for v in test_result.values()]) == 43
    assert np.array_equal(
        test_result[(19,21)], [np.array(( 0.03322,  0.23993,  0.12885))])
    assert np.array_equal(
        test_result[(1,2)]  , [np.array(( 0.27009,  0.03611, -0.16173))])
    assert np.array_equal(
        test_result[(29,24)], [np.array(( 0.05682,  0.09902, -0.04100))])
    assert np.array_equal(
        test_result[(29,30)], [np.array((-0.06080, -0.04975, -0.19957))])
    assert np.array_equal(
        test_result[(13,12)], [
        np.array((-0.18853, -0.19153,  0.12901)),
        np.array((-0.36498, -1.01171,  0.10899)),
        np.array((-0.63112, -1.67971,  0.29640)),
    ])

def test_get_bond_dipole():
    """as name."""
    test_stru = sp.get_structure(f"{STRU_DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J": (0, 1)})
    test_region = create_region_from_selection_pattern(test_stru, "resi 101+254")
    target_bond = (test_stru.ligands[0].find_atom_name("CAE"), test_stru.ligands[0].find_atom_name("H2"))
    test_ele_stru = ElectronicStructure(
        energy_0=0.0,
        geometry=test_region,
        mo=f"{DATA_DIR}../../analysis/data/KE_07_R7_2_S_101_254.chk",
        mo_parser=None,
        source="gaussian16",
    )
    cluster_job_config = {
        "cluster": Accre(),
        "res_keywords": {
            "account": "yang_lab_csb",
            "partition": "debug",
            'walltime': '30:00',
        }
    }
    result = mi.get_bond_dipole(
        test_ele_stru,
        target_bond[0],
        target_bond[1],
        work_dir=WORK_DIR,
        keep_in_file=False,
        cluster_job_config=cluster_job_config,
        job_check_period=3,
    )

    assert result == (0.39995, np.array(0.33744, -0.21187, -0.03473))
