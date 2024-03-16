"""Testing enzy_htp.analysis.stability.py

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-03-12"""
import os
import numpy as np

from enzy_htp.analysis import ddg_fold_of_mutants
from enzy_htp import PDBParser
from enzy_htp.mutation import assign_mutant
from enzy_htp.core.clusters.accre import Accre

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()

def test_ddg_fold_of_mutants():
    """as name.
    use data from EnzyHTP 1.0 as answer"""
    os.environ["ROSETTA3"] = "/data/yang_lab/Common_Software/Rosetta3.9/main"
    test_stru = sp.get_structure(f"{STRU_DATA_DIR}KE_07_R7_2_S.pdb")
    mutant_space = assign_mutant(
        test_stru,
        "{E24V,K162I,R163L},{S29K,E24V,K162L,R163L}"
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "partition" : "production",
            "account" : "yang_lab",
        }
    }
    result = ddg_fold_of_mutants(
        test_stru,
        mutant_space,
        num_iter = 3,
        cluster_job_config = cluster_job_config,
        relax_cluster_job_config = cluster_job_config,
        work_dir=WORK_DIR,
        action_on_wt_at_start=False, # Turn off relax temp for testing
    )

    mut_1 = tuple(mutant_space[0])
    mut_2 = tuple(mutant_space[1])
    assert np.isclose(result[mut_1], -0.142, atol=1)
    assert np.isclose(result[mut_2], -1.2, atol=1)

