"""Testing enzy_htp.geometry.sampling.py
Author
Date
"""

import glob
import pytest
import os

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.structure import structure_constraint
from enzy_htp.geometry import md_simulation, equi_md_sampling
from enzy_htp import interface
from enzy_htp import PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()
amber_interface = interface.amber


# TODO: finish these tests while finished Amber interface
@pytest.mark.accre
@pytest.mark.long
def test_md_simulation_amber_lv1():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine
    level 1:
    - no replica
    - no constrain"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }
    step_1  = amber_interface.build_md_step(
        name="min",
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",)

    step_2 = amber_interface.build_md_step(
        name="equi_npt",
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,)

    step_3 = amber_interface.build_md_step(
        name="prod_npt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        restart=True,
        if_report=True,
        record_period=0.0005,) # TODO address keep_in_file

    md_result = md_simulation(stru=test_stru,
                  param_method=test_param_method,
                  steps=[step_1, step_2, step_3],
                  parallel_runs=1,
                  job_check_period=10)

    # clean up
    fs.clean_temp_file_n_dir([

    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob("MD/rep_0/*out")
    + glob.glob("MD/rep_0/*nc")
    + glob.glob("MD/rep_0/*rst")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))

@pytest.mark.accre
def test_md_simulation_amber_lv2():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine
    level 1:
    - no replica
    - backbone constrain"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{DATA_DIR}/ncaa_lib_empty",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "period" : 60,
        "res_keywords" : {"account" : "csb_gpu_acc"}
    }
    constrain = structure_constraint.build_from_preset(test_stru, "freeze_backbone")
    step_1  = amber_interface.build_md_step(
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="GPU",
        constrain=constrain,)

    step_2 = amber_interface.build_md_step(
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="GPU",
        temperature=300,
        constrain=constrain,)

    step_3 = amber_interface.build_md_step(
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="GPU",
        if_report=True,
        temperature=300,
        record_period=0.0005,)

    md_simulation(stru=test_stru,
                  param_method=test_param_method,
                  steps=[step_1, step_2, step_3],
                  parallel_runs=1)

@pytest.mark.accre
def test_md_simulation_amber_3_repeat():
    """Test running a 3-replica MD.
    Using Amber & Accre as an example engine"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_param_method = amber_interface.TODO
    step_1 = amber_interface.TODO
    step_2 = amber_interface.TODO
    step_3 = amber_interface.TODO

    md_simulation(stru=test_stru,
                  param_method=test_param_method,
                  steps=[step_1, step_2, step_3],
                  parallel_runs=3)

@pytest.mark.accre
def test_equi_md_sampling_lv1():
    """test for equi_md_sampling
    level 1:"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)}) # this should be explicitly set
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{DATA_DIR}/ncaa_lib_empty",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "period" : 600,
        "res_keywords" : {"account" : "csb_gpu_acc"}
    }
    equi_md_sampling(stru = test_stru,
                     param_method = test_param_method,
                     cluster_job_config = cluster_job_config,)


@pytest.mark.accre
def test_equi_md_sampling_lv2():
    """test for equi_md_sampling
    level 2:
    - geom constrain"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{DATA_DIR}/ncaa_lib_empty",
    )
    test_prod_constrain = None
    cluster_job_config = {
        "cluster" : Accre(),
        "period" : 600,
        "res_keywords" : {"account" : "csb_gpu_acc"}
    }
    equi_md_sampling(stru = test_stru,
                     param_method = test_param_method,
                     prod_constrain = test_prod_constrain,
                     cluster_job_config = cluster_job_config,)
