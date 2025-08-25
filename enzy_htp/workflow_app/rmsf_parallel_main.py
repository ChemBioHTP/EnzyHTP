"""This script calculates RMSF for trajs in parallel.
NOTE: PyMol is not thread safe so we can either:
1. limit threads in each subprocesses to 1 (as done here)
or 
2. re-import everything in unit_task and limit threads in processes to 1

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-11
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import multiprocessing
import logging
from pathlib import Path

from enzy_htp import interface as eh_interface
from enzy_htp.analysis import rmsf
from enzy_htp.core.general import save_obj
from enzy_htp.core.logger import _LOGGER as eh_logger
from enzy_htp import config as eh_config
import enzy_htp.core.file_system as fs

def calculate_rmsf_from_traj(traj_file: str, top_file: str, ref_pdb: str):
    """calculate RMSF for each residues from a trajectorie"""
    structure_ensemble = eh_interface.amber.load_traj(
        prmtop_path=top_file,
        traj_path=traj_file,
        ref_pdb=ref_pdb,
    )
    result = rmsf(stru_esm=structure_ensemble, by_residue=True)
    fs.safe_rm(structure_ensemble.topology_source_file)
    return list(result.values())

def unit_task(args):
    """unit task of the parallel run"""
    idx, top_file, traj_file, ref_pdb, group, scratch = args

    eh_config.system.SCRATCH_DIR = scratch # avoid non-lock files to conflict

    result = {}
    result_file = group / "rmsf.pickle"
    print(f"working on {idx}")
    if traj_file.exists():
        if not result_file.exists():
            result[int(idx)] = calculate_rmsf_from_traj(
                top_file = str(top_file), traj_file = str(traj_file), ref_pdb = str(ref_pdb))
            print(f"finished the calculation for {idx}")
            save_obj(result, str(result_file))
        else:
            print(f"{result_file} already exists")
    else:
        print(f"MD failed in {idx}")

def create_rmsf_dataset():
    """calculate RMSF for each residues from trajectories in the dataset"""

    task_args = []
    eh_logger.setLevel(logging.WARNING)

    for master_dir in [
        "data/training_set/training_set_md/",
        "data/test_set/test_set_md/",
    ]:
        groups = Path(master_dir).glob("group_*")
        for group in groups:
            top_file = next(group.glob("traj_50/*prmtop"))
            traj_file = group / "traj_50/MD/prod.mdcrd"
            ref_pdb = next(group.glob("traj_50/*_aH.pdb"))
            idx = ref_pdb.stem.removesuffix("_rmW_rmL_rmH_aH")
            scratch = group  / "scratch"
            fs.safe_mkdir(str(scratch))
            task_args.append((idx, top_file, traj_file, ref_pdb, group, scratch))

    with multiprocessing.Pool(processes=40) as pool:
        pool.map(unit_task, task_args)

def main():
    create_rmsf_dataset()

if __name__ == "__main__":
    main()
