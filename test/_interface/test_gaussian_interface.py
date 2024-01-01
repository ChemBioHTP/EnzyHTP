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
from enzy_htp.structure.structure_region import create_region_from_selection_pattern
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
    assert md_step.capping_method == "res_ter_cap"

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

def test_get_method_keyword_from_name():
    """as name"""
    test_name = "b3lyp-d3"
    answer_kw = "b3lyp em=gd3"
    test_kw = gi.get_method_keyword_from_name(test_name)
    assert answer_kw == test_kw

    test_name = "b3lyp-d3bj"
    answer_kw = "b3lyp em=gd3bj"
    test_kw = gi.get_method_keyword_from_name(test_name)
    assert answer_kw == test_kw

    test_name = "b3lyp-D3(BJ)"
    answer_kw = "b3lyp em=gd3bj"
    test_kw = gi.get_method_keyword_from_name(test_name)
    assert answer_kw == test_kw

    test_name = "pbe0-D3(BJ)"
    answer_kw = "pbe1pbe em=gd3bj"
    test_kw = gi.get_method_keyword_from_name(test_name)
    assert answer_kw == test_kw

def test_lot_to_keyword():
    """as name"""
    test_lot = QMLevelofTheory(
        basis_set="6-31G(d)",
        method="HF",
        solvent="water",
        solv_method="SMD",
    )
    answer_kw = ("HF 6-31G(d) scrf=(SMD, solvent=water)", [], [])
    test_kw = gi.lot_to_keyword(test_lot)
    assert answer_kw == test_kw

def test_get_geom_lines():
    """as name.
    answer confirmed using GaussView manually"""
    test_atoms = sp.get_structure(f"{DATA_DIR}H5J.pdb").atoms
    answer_lines = [
        "O         26.31000000     -2.47500000    -48.71100000",
        "N         26.99300000     -2.19100000    -49.70300000",
        "O         28.21800000     -1.97000000    -49.63900000",
        "C         26.34700000     -2.11100000    -50.99500000",
        "C         24.96800000     -2.35600000    -51.09900000",
        "H         24.35300000     -2.60600000    -50.24700000",
        "C         27.13100000     -1.78700000    -52.12100000",
        "H         28.18900000     -1.61100000    -51.99600000",
        "C         26.56500000     -1.69700000    -53.38000000",
        "H         27.16700000     -1.44800000    -54.24200000",
        "C         25.17500000     -1.93700000    -53.52700000",
        "C         24.40200000     -2.26500000    -52.36300000",
        "C         23.03800000     -2.45000000    -52.83000000",
        "H         22.29200000     -2.70700000    -52.09200000",
        "N         22.84300000     -2.29500000    -54.05600000",
        "O         24.50900000     -1.89000000    -54.63400000",
    ]
    geom_lines = gi.get_geom_lines(test_atoms, [])

    assert geom_lines == answer_lines

def test_get_geom_lines_w_cons(): # TODO finish this when needed
    """as name.
    answer confirmed using GaussView manually"""
    test_atoms = sp.get_structure(f"{DATA_DIR}H5J.pdb").atoms
    answer_lines = [
    ]
    geom_lines = gi.get_geom_lines(test_atoms, ["TODO"])
    assert False

def test_make_mol_spec():
    """as name"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_region = create_region_from_selection_pattern(
        "br. (resi 254 around 5)", test_stru
    )
    constraints = []
    answer = [

    ]
    test_result = gi._make_mol_spec(
        test_stru,
        test_stru_region,
        constraints,
    )
    print(*test_result, sep="\n")
    # assert test_result == answer
