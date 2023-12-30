"""Testing enzy_htp.geometry.sampling.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-25
"""

import glob
import pytest
import os

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp.geometry import md_simulation, equi_md_sampling
from enzy_htp import interface
from enzy_htp import PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()
amber_interface = interface.amber


@pytest.mark.accre
@pytest.mark.long
def test_md_simulation_amber_lv1():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine
    level 1:
    - no replica
    - no constrain
    took around 2 min. use the internal checking in translate() to
    make sure each step is successfully finished."""
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
        record_period=0.0005,)

    md_result = md_simulation(stru=test_stru,
                  param_method=test_param_method,
                  steps=[step_1, step_2, step_3],
                  parallel_runs=1,
                  job_check_period=10)

    # TODO may be also check MD traj is reasonable?

    # clean up
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob("MD/rep_0/*out")
    + glob.glob("MD/rep_0/*nc")
    + glob.glob("MD/rep_0/*rst")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))

@pytest.mark.accre
@pytest.mark.long
def test_md_simulation_amber_lv2():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine
    level 1:
    - no replica
    - backbone constrain"""
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

    constrain = [stru_cons.create_backbone_freeze(test_stru)]

    step_1  = amber_interface.build_md_step(
        name="min",
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=constrain,)

    step_2 = amber_interface.build_md_step(
        name="equi_npt",
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        constrain=constrain,)

    step_3 = amber_interface.build_md_step(
        name="prod_npt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        restart=True,
        if_report=True,
        record_period=0.0005,
        constrain=constrain,)

    md_result = md_simulation(
        stru=test_stru,
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
@pytest.mark.long
def test_md_simulation_amber_lv3():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine
    level 1:
    - no replica
    - backbone constrain
    - distance and angle constrain"""
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

    constrain = [stru_cons.create_backbone_freeze(test_stru),
                 stru_cons.create_distance_constraint(
                    "B.254.H2", "A.101.OE2", 2.4, test_stru),
                 stru_cons.create_angle_constraint(
                    "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru),]

    step_1  = amber_interface.build_md_step(
        name="min",
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=constrain,)

    step_2 = amber_interface.build_md_step(
        name="equi_npt",
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        constrain=constrain,)

    step_3 = amber_interface.build_md_step(
        name="prod_npt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        restart=True,
        if_report=True,
        record_period=0.0005,
        constrain=constrain,)

    md_result = md_simulation(
        stru=test_stru,
        param_method=test_param_method,
        steps=[step_1, step_2, step_3],
        parallel_runs=1,
        job_check_period=10)

    assert len(glob.glob("MD/rep_0/*in")) == 0

    # clean up
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob("MD/rep_0/*out")
    + glob.glob("MD/rep_0/*nc")
    + glob.glob("MD/rep_0/*rst")
    + glob.glob("MD/rep_0/*rs")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))


@pytest.mark.accre
@pytest.mark.long
def test_md_simulation_amber_3_repeat():
    """Test running a 3-replica MD.
    Using Amber & Accre as an example engine"""
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

    constrain = [stru_cons.create_backbone_freeze(test_stru)]

    step_1  = amber_interface.build_md_step(
        name="min",
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=constrain,)

    step_2 = amber_interface.build_md_step(
        name="equi_npt",
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        constrain=constrain,)

    step_3 = amber_interface.build_md_step(
        name="prod_npt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        restart=True,
        if_report=True,
        record_period=0.0005,
        constrain=constrain,)

    md_result = md_simulation(
        stru=test_stru,
        param_method=test_param_method,
        steps=[step_1, step_2, step_3],
        parallel_runs=3,
        job_check_period=10)

    # clean up
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob("MD/rep_0/*out")
    + glob.glob("MD/rep_0/*nc")
    + glob.glob("MD/rep_0/*rst")
    + glob.glob("MD/rep_1/*out")
    + glob.glob("MD/rep_1/*nc")
    + glob.glob("MD/rep_1/*rst")
    + glob.glob("MD/rep_2/*out")
    + glob.glob("MD/rep_2/*nc")
    + glob.glob("MD/rep_2/*rst")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))

@pytest.mark.accre
@pytest.mark.long
@pytest.mark.temp
def test_equi_md_sampling_lv1():
    """test for equi_md_sampling
    level 1: no constraint"""
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
    md_result = equi_md_sampling(
        stru = test_stru,
        param_method = test_param_method,
        cluster_job_config = cluster_job_config,
        job_check_period=10,
        # shorter sim for test
        prod_time=0.5,
        record_period=0.05,
        work_dir=f"{WORK_DIR}MD/")

    # clean up
    # fs.safe_rmdir(f"{WORK_DIR}MD/")
    # fs.clean_temp_file_n_dir([
    # ] + glob.glob("slurm-*.out")
    # + glob.glob("scratch/amber_parameterizer/*")
    # + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))


@pytest.mark.accre
def test_equi_md_sampling_lv2():
    """test for equi_md_sampling
    level 2:
    - geom constrain"""
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
    constrain = [stru_cons.create_distance_constraint(
                    "B.254.H2", "A.101.OE2", 2.4, test_stru),
                 stru_cons.create_angle_constraint(
                    "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru),]
    md_result = equi_md_sampling(
        stru = test_stru,
        param_method = test_param_method,
        cluster_job_config = cluster_job_config,
        job_check_period=10,
        prod_constrain=constrain,
        # shorter sim for test
        prod_time=0.5,
        record_period=0.05,
        work_dir=f"{WORK_DIR}MD/")

    # clean up
    fs.safe_rmdir(f"{WORK_DIR}MD/")
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))

def test_equi_md_sampling_wrong_cons():
    """test for equi_md_sampling that uses wrong constraints"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_2 = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S_mut.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }
    constrain = [stru_cons.create_distance_constraint(
                    "B.254.H2", "A.101.OE2", 2.4, test_stru_2),
                 stru_cons.create_angle_constraint(
                    "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru_2),]

    with pytest.raises(ValueError) as e:
        md_result = equi_md_sampling(
            stru = test_stru,
            param_method = test_param_method,
            cluster_job_config = cluster_job_config,
            job_check_period=10,
            prod_constrain=constrain,
            # shorter sim for test
            prod_time=0.5,
            record_period=0.05,
            work_dir=f"{WORK_DIR}MD/")
