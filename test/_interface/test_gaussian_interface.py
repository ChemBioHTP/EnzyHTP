"""Testing enzy_htp._interface.gaussian_interface.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-12-29
"""
import pytest
import os

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.chemical.level_of_theory import QMLevelofTheory, MMLevelofTheory
from enzy_htp._config.armer_config import ARMerConfig
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp import PDBParser
from enzy_htp import interface
from enzy_htp._interface.gaussian_interface import GaussianSinglePointEngine

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()
gi = interface.gaussian

@pytest.mark.accre
def test_single_point_make_job_lv1():
    """Test GaussianSinglePointEngine().make_job()
    Using Accre as an example cluster
    level 1:
    - full structure (only 1 ligand)"""
    test_stru = sp.get_structure(f"{DATA_DIR}H5J.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_method = QMLevelofTheory(
        basis_set="6-31G(d)",
        method="HF",
        solvent="water",
        solv_method="SMD",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "yang_lab_csb",
                         "partition" : "production"}
    }
    
    qm_engine = gi.build_single_point_engine(
        region=[],
        method=test_method,
        keep_geom=True,
        cluster_job_config=cluster_job_config,
    )
    test_job, test_egg = qm_engine.make_job(test_stru)

def test_build_single_point_engine_default():
    """as said in the name. Assert several default values
    as samples."""
    gi = interface.gaussian
    md_step: GaussianSinglePointEngine = gi.build_single_point_engine()
    assert md_step.method == QMLevelofTheory(
        basis_set = "def2-svp",
        method = "pbe0",
        solvent = None,
        solv_method = None,
    )
    assert md_step.region == None
    assert md_step.cluster_job_config["res_keywords"] == ARMerConfig.QM_SPE_CPU_RES
    assert md_step.keep_geom == True
    assert md_step.capping_method == "h_cap"

def test_build_single_point_engine_res_keywords():
    """as said in the name. Assert several default values
    as samples."""
    gi = interface.gaussian
    md_step: GaussianSinglePointEngine = gi.build_single_point_engine(
        cluster_job_config={
            "cluster" : None,
            "res_keywords" : {"partition" : "production",
                                "account" : "yang_lab",}
        })
    assert md_step.cluster_job_config["cluster"] is None
    assert md_step.cluster_job_config["res_keywords"]["core_type"] == "cpu"
    assert md_step.cluster_job_config["res_keywords"]["partition"] == "production"
    assert md_step.cluster_job_config["res_keywords"]["nodes"] == "1"
    assert md_step.cluster_job_config["res_keywords"]["node_cores"] ==  "8"
    assert md_step.cluster_job_config["res_keywords"]["job_name"] ==  "QM_SPE_EnzyHTP"
    assert md_step.cluster_job_config["res_keywords"]["mem_per_core"] ==  "3G"
    assert md_step.cluster_job_config["res_keywords"]["walltime"] ==  "3-00:00:00"
    assert md_step.cluster_job_config["res_keywords"]["account"] ==  "yang_lab"
