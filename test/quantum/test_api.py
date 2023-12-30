"""Testing enzy_htp.quantum.api.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-12-29
"""

import pytest
import os

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.chemical.level_of_theory import QMLevelofTheory, MMLevelofTheory
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp import interface
from enzy_htp import PDBParser
from enzy_htp.quantum import single_point

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
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
    test_method = QMLevelofTheory(
        basis_set="def2-SVP",
        method="PBE0",
        solvent="water",
        solv_method="SMD",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "yang_lab_csb",
                         "partition" : "production"}
    }
    
    qm_result = single_point(
        stru=test_stru,
        engine="gaussian",
        method=test_method,
        cluster_job_config=cluster_job_config,
        job_check_period=10,
    )

@pytest.mark.accre
def test_single_point_gaussian_lv2():
    """Test running a QM single point.
    Using Gaussian & Accre as an example engine
    level 2:
    - structure region
    - single snapshot"""
    test_stru = sp.get_structure(f"{DATA_DIR}H5J.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_method = QMLevelofTheory(
        basis_set="def2-SVP",
        method="PBE0",
        solvent="water",
        solv_method="SMD",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "yang_lab_csb",
                         "partition" : "production"}
    }
    
    qm_result = single_point(
        stru=test_stru,
        engine="gaussian",
        method=test_method,
        cluster_job_config=cluster_job_config,
        job_check_period=10,
        )

@pytest.mark.accre
def test_single_point_gaussian_lv3():
    """Test running a QM single point.
    Using Gaussian & Accre as an example engine
    level 3:
    - 1 structure region
    - 100 snapshots"""
    test_stru = sp.get_structure(f"{DATA_DIR}H5J.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_method = QMLevelofTheory(
        basis_set="def2-SVP",
        method="PBE0",
        solvent="water",
        solv_method="SMD",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "yang_lab_csb",
                         "partition" : "production"}
    }
    
    qm_result = single_point(
        stru=test_stru,
        engine="gaussian",
        method=test_method,
        cluster_job_config=cluster_job_config,
        job_check_period=10,
        )

@pytest.mark.accre
def test_single_point_gaussian_lv4():
    """Test running a QM single point.
    Using Gaussian & Accre as an example engine
    level 4:
    - 2 structure region QM/MM
    - 100 snapshots"""
    test_stru = sp.get_structure(f"{DATA_DIR}H5J.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_method = QMLevelofTheory(
        basis_set="def2-SVP",
        method="PBE0",
        solvent="water",
        solv_method="SMD",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "yang_lab_csb",
                         "partition" : "production"}
    }
    
    qm_result = single_point(
        stru=test_stru,
        engine="gaussian",
        method=test_method,
        cluster_job_config=cluster_job_config,
        job_check_period=10,
        )



