"""Testing enzy_htp.analysis.clustering.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2025-07-08
"""
import os
import numpy as np

from enzy_htp.core.clusters.accre import Accre
from enzy_htp.analysis import rep_stru_clustering
from enzy_htp import PDBParser, interface

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
INT_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../_interface/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()

def test_rep_stru_clustering():
    """Test running bond_dipole calculation
    Using kmeans_rmsd as an example set up"""
    stru_esm = interface.amber.load_traj(
        prmtop_path=f"{INT_DATA_DIR}/mmpbsa_test_sol.prmtop",
        traj_path=f"{INT_DATA_DIR}/mmpbsa_test_sol_10f.nc",
        ref_pdb=f"{INT_DATA_DIR}/mmpbsa_test_sol.pdb"
    )
    ligand = "resi 290"
    answer = [
        (0.3, np.array([3, 8, 9])),
        (0.2, np.array([1, 4])),
        (0.5, np.array([0, 2, 5, 6, 7])),
    ]
    
    result = rep_stru_clustering(
        stru_esm, ligand,
        method="kmeans_rmsd",
        n_clusters = 3, random_state = 42,
    )

    for r, a in zip(result.values(), answer):
        assert r[0] == a[0]
        assert r[1].all() == a[1].all()
