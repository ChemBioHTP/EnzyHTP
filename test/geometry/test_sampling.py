"""Testing enzy_htp.geometry.sampling.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-25
"""

import pytest
import os

from enzy_htp.geometry import md_simulation
from enzy_htp import interface
from enzy_htp import PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
sp = PDBParser()
amber_interface = interface.amber


# TODO: finish these tests while finished Amber interface EOD
@pytest.mark.accre
def test_md_simulation_amber_no_repeat():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_parallel_method = amber_interface.TODO
    step_1 = amber_interface.TODO
    step_2 = amber_interface.TODO
    step_3 = amber_interface.TODO

    md_simulation(stru=test_stru,
                  parallel_method=test_parallel_method,
                  steps=[step_1, step_2, step_3],
                  parallel_runs=1)

@pytest.mark.accre
def test_md_simulation_amber_no_repeat():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_parallel_method = amber_interface.TODO
    step_1 = amber_interface.TODO
    step_2 = amber_interface.TODO
    step_3 = amber_interface.TODO

    md_simulation(stru=test_stru,
                  parallel_method=test_parallel_method,
                  steps=[step_1, step_2, step_3],
                  parallel_runs=3)
