"""Testing enzy_htp.analysis.binding.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-11-07
"""
import os
import numpy as np

from enzy_htp.core.clusters.accre import Accre
from enzy_htp.analysis import binding_energy
from enzy_htp import PDBParser, interface

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
INT_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../_interface/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()

def test_binding_energy():
    """Test running bond_dipole calculation
    Using LMO-Multiwfn & Accre as an example set up
    
    Note: need 4 CPUs to run this."""
    stru_esm = interface.amber.load_traj(
        prmtop_path=f"{INT_DATA_DIR}/mmpbsa_test_sol.prmtop",
        traj_path=f"{INT_DATA_DIR}/mmpbsa_test_sol_10f.nc",
        ref_pdb=f"{INT_DATA_DIR}/mmpbsa_test_sol.pdb"
    )
    ligand = "resi 290"
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "partition" : "production",
            "account" : "yang_lab",}
    }
    
    result = binding_energy(
        stru_esm, ligand,
        cluster_job_config=cluster_job_config,
    )

    assert np.isclose(result, -1.9344, atol=1)
