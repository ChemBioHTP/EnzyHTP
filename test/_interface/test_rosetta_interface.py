"""Testing enzy_htp._interface.rosetta_interface.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-02-16
"""
from pathlib import Path
import pickle
import re
import numpy as np
import pytest
import os

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp._config.armer_config import ARMerConfig
from enzy_htp import PDBParser
from enzy_htp import interface, config
from enzy_htp.core.job_manager import ClusterJob
from enzy_htp.mutation_class.mutation import Mutation
from enzy_htp.structure.structure import Structure

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()
ri = interface.rosetta



def test_relax_w_stru():
    """as said in the name. make sure relax works Structure() input."""
    os.environ["ROSETTA3"] = "/data/yang_lab/Common_Software/Rosetta3.9/main"
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    scores = ri.relax(
        test_stru,
        nstruct=3,
        constrain_relax_to_start_coords=False,
        coord_constrain_sidechains=False,
        ramp_constraints=False,
        use_cartesian=True,
        scorefxn="ref2015_cart",
        cluster_job_config={
            "cluster" : Accre(),
            "res_keywords" : {
                "partition" : "production",
                "node_cores" : "24",
                "account" : "yang_lab",
                "walltime" : "10:00:00",
            }
        },
        delete_crash=False,
    )
    assert len(scores) == 3
    for pdb_file in scores["description"]:
        assert Path(pdb_file).exists()

def test_get_structure_from_score():
    """as said in the name. use a hand made score_df in a pickle file.
    just make sure the Structure is returned for now"""
    with open(f"{DATA_DIR}test_rosetta_relax_score/score_df.pickle", "rb") as f:
        test_score_df = pickle.load(f)
    result = ri.get_structure_from_score(test_score_df, base_dir=".", clean_up_pdb=False)
    assert isinstance(result, Structure)

def test_get_ddg_fold():
    """as said in the name."""
    test_ddg_file = f"{DATA_DIR}test.ddg"
    ddg_engine = ri.build_cartesian_ddg_engine()
    result = ddg_engine.get_ddg_fold(test_ddg_file)
    assert np.isclose(result, 12.069600000000037, atol=1e-6)

def test_cartesian_ddg_action_on_wt():
    """as said in the name.
    just make sure the Structure is returned for now"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    ddg_engine = ri.build_cartesian_ddg_engine(
        relax_cluster_job_config={
            "cluster" : Accre(),
            "res_keywords" : {
                "partition" : "production",
                "account" : "yang_lab",
            }
        }
    )
    new_stru = ddg_engine.action_on_wt(test_stru)
    PDBParser().save_structure(f"{DATA_DIR}KE_07_R7_2_S_cart_relax.pdb", new_stru)
    # assert isinstance(new_stru, Structure)

def test_cartesian_ddg_make_job(helpers):
    """as said in the name.
    just make sure no error occurs for now"""
    os.environ["ROSETTA3"] = "/data/yang_lab/Common_Software/Rosetta3.9/main"
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S_cart_relax.pdb")
    ddg_engine = ri.build_cartesian_ddg_engine(
        cluster_job_config={
            "cluster" : Accre(),
            "res_keywords" : {
                "partition" : "production",
                "account" : "yang_lab",
            }
        },
        work_dir=WORK_DIR
    )
    mutant = [
        Mutation(orig='ARG', target='TRP', chain_id='A', res_idx=154),
        Mutation(orig='HIS', target='ALA', chain_id='A', res_idx=201),
        ]
    job, result_egg = ddg_engine.make_job(test_stru, mutant)

    assert helpers.equiv_files(f"{WORK_DIR}/mutations.txt", f"{DATA_DIR}/mutations_answer.txt")
    # TODO assert on job.sub_script_str
    
    fs.safe_rm(f"{WORK_DIR}/mutations.txt")
    fs.safe_rm(f"{WORK_DIR}/cart_ddg_temp.pdb")

def test_check_dgg_fold_error():
    """as said in the name.
    use a pseudo ClusterJob. Only testing for ddg file existence
    just make sure no error occurs for now"""
    os.environ["ROSETTA3"] = "/data/yang_lab/Common_Software/Rosetta3.9/main"
    test_ddg_file = f'{DATA_DIR}mutations.ddg'
    test_job = ClusterJob(Accre(), "")
    test_job.job_cluster_log = f"{DATA_DIR}empty_file"
    ddg_engine = ri.build_cartesian_ddg_engine(
        cluster_job_config={
            "cluster" : Accre(),
            "res_keywords" : {
                "partition" : "production",
                "account" : "yang_lab",
            }
        },
        work_dir=WORK_DIR
    )
    ddg_engine.check_dgg_fold_error(
        test_ddg_file, test_job
        )
