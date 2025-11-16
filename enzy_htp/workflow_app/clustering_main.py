import logging

from enzy_htp.analysis import rep_stru_clustering
from enzy_htp import interface
from enzy_htp import PDBParser
import enzy_htp.core.file_system as fs
from enzy_htp import _LOGGER

MD_DIR = "MD/mutant_0/MD"
sp = PDBParser()

stru_esm = interface.amber.load_traj(
    prmtop_path=f"{MD_DIR}/rmwat.prmtop",
    traj_path=f"{MD_DIR}/prod_npt_merged_rmwat.nc",
)
ligand = "resn LIG"

result = rep_stru_clustering(
    stru_esm, ligand,
    method="kmeans_rmsd",
    n_clusters = 8,
)

result_dir = "clustering/WT"
fs.safe_mkdir(result_dir)

for i, (k, v) in enumerate(result.items()):
    weight, members = v
    sp.save_structure(f"{result_dir}/rep_{i}_{weight}.pdb", k)
