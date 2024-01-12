"""Testing enzy_htp.quantum.api.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-12-29
"""

import numpy as np
import pytest
import os
from enzy_htp._interface.amber_interface import AmberMDCRDParser

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.chemical.level_of_theory import QMLevelOfTheory, MMLevelOfTheory
from enzy_htp.core.general import get_itself
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp import interface
from enzy_htp import PDBParser
from enzy_htp.quantum import single_point
from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.structure.structure_region.api import create_region_from_selection_pattern

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../structure/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()
gi = interface.gaussian


@pytest.mark.accre
def test_single_point_gaussian_lv1():
    """Test running a QM single point.
    Using Gaussian & Accre as an example engine
    level 1:
    - full structure (only 1 ligand)
    - single snapshot"""
    test_stru = sp.get_structure(f"{DATA_DIR}H5J.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_method = QMLevelOfTheory(
        basis_set="3-21G",
        method="HF",
        solvent="water",
        solv_method="SMD",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "account" : "yang_lab_csb",
            "partition" : "debug",
            'walltime' : '30:00',
        }
    }

    qm_result = single_point(
        stru=test_stru,
        engine="gaussian",
        method=test_method,
        cluster_job_config=cluster_job_config,
        job_check_period=10,
        work_dir=f"{WORK_DIR}/QM_SPE/"
    )[0]

    assert np.isclose(qm_result.energy_0, -597.293275805, atol=1e-3)
    assert not os.path.exists(f"{WORK_DIR}/QM_SPE/gaussian_spe.gjf")

    fs.safe_rm(f"{WORK_DIR}/QM_SPE/gaussian_spe.out")
    fs.safe_rm(f"{WORK_DIR}/QM_SPE/gaussian_spe.chk")
    fs.safe_rmdir(f"{WORK_DIR}/QM_SPE/", True)

@pytest.mark.accre
def test_single_point_gaussian_lv2():
    """Test running a QM single point.
    Using Gaussian & Accre as an example engine
    level 2:
    - structure region
    - single snapshot"""
    test_stru = sp.get_structure(f"{STRU_DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_method = QMLevelOfTheory(
        basis_set="3-21G",
        method="HF",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "account" : "yang_lab_csb",
            "partition" : "debug", # cant use debug that accre only have a node that dont have Gaussian
            'walltime' : '30:00',
        }
    }

    qm_result = single_point(
        stru=test_stru,
        engine="gaussian",
        method=test_method,
        regions=["br. (resi 254 around 2)"],
        cluster_job_config=cluster_job_config,
        job_check_period=10,
        work_dir=f"{WORK_DIR}/QM_SPE/",
        )[0]

    assert np.isclose(qm_result.energy_0, -1668.68154951, atol=1e-3)
    assert not os.path.exists(f"{WORK_DIR}/QM_SPE/gaussian_spe.gjf")

    fs.safe_rm(f"{WORK_DIR}/QM_SPE/gaussian_spe.out")
    fs.safe_rm(f"{WORK_DIR}/QM_SPE/gaussian_spe.chk")
    fs.safe_rmdir(f"{WORK_DIR}/QM_SPE/", True)

@pytest.mark.accre
@pytest.mark.temp
def test_single_point_gaussian_lv3():
    """Test running a QM single point.
    Using Gaussian & Accre as an example engine
    level 3:
    - 1 structure region
    - 100 snapshots"""
    test_stru = sp.get_structure(f"{STRU_DATA_DIR}KE_07_R7_2_S.pdb")
    test_traj = f"{STRU_DATA_DIR}KE_07_R7_2_S_10f.mdcrd"
    test_prmtop = f"{STRU_DATA_DIR}KE_07_R7_2_S_10f.prmtop"
    test_esm = StructureEnsemble(
        topology=test_stru,
        top_parser=get_itself,
        coordinate_list=test_traj,
        coord_parser=AmberMDCRDParser(test_prmtop).get_coordinates,
    )
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_method = QMLevelOfTheory(
        basis_set="3-21G",
        method="HF",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "account" : "yang_lab_csb",
            "partition" : "debug",
            'walltime' : '30:00',
        }
    }

    qm_results = single_point(
        stru=test_esm,
        engine="gaussian",
        method=test_method,
        regions=["br. (resi 254 around 2)"],
        cluster_job_config=cluster_job_config,
        job_check_period=10,
        job_array_size=5,
        work_dir=f"{WORK_DIR}/QM_SPE/",
        )

    for qm_result in qm_results:
        assert np.isclose(qm_result.energy_0, -1668, atol=5)
    fs.safe_rmdir(f"{WORK_DIR}/QM_SPE/")

@pytest.mark.accre
def test_single_point_gaussian_lv4():
    """Test running a QM single point.
    Using Gaussian & Accre as an example engine
    level 4:
    - 2 structure region QM/MM
    - 100 snapshots"""
    assert False
    test_stru = sp.get_structure(f"{DATA_DIR}H5J.pdb")
    test_stru.assign_ncaa_chargespin({"H5J": (0, 1)})
    test_method = QMLevelOfTheory(
        basis_set="3-21G",
        method="HF",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "account" : "yang_lab_csb",
            "partition" : "debug",
            'walltime' : '30:00',
        }
    }


    qm_result = single_point(
        stru=test_stru,
        engine="gaussian",
        method=test_method,
        cluster_job_config=cluster_job_config,
        job_check_period=10,
    )
