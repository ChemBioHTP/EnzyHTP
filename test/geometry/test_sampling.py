"""Testing enzy_htp.geometry.sampling.py
Author
Date
"""

import pytest
import os

from enzy_htp.core.clusters.accre import Accre
from enzy_htp.geometry import md_simulation, equi_md_sampling
from enzy_htp import interface
from enzy_htp import PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
sp = PDBParser()
amber_interface = interface.amber


# TODO: finish these tests while finished Amber interface
@pytest.mark.accre
def test_md_simulation_amber_no_repeat():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{DATA_DIR}/ncaa_lib_empty",
    )
    step_1 = amber_interface.build_md_step()
    step_2 = amber_interface.build_md_step()
    step_3 = amber_interface.build_md_step()

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
        "res_setting" : {"account" : "csb_gpu_acc"}
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
        "res_setting" : {"account" : "csb_gpu_acc"}
    }
    equi_md_sampling(stru = test_stru,
                     param_method = test_param_method,
                     prod_constrain = test_prod_constrain,
                     cluster_job_config = cluster_job_config,)
