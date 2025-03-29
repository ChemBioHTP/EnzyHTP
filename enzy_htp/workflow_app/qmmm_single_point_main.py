"""This is MVP of running a Gaussian/Amber based QMMM single point for every frame in a trajectory.
TODO convert it into an engine of the science API"""

import csv
import glob
import os
from typing import Dict, List

from enzy_htp.core.clusters.accre import Accre
from enzy_htp.core.job_manager import ClusterJob
import enzy_htp.core.file_system as fs
from enzy_htp.core.general import save_obj, load_obj
from enzy_htp.chemical import QMLevelOfTheory
from enzy_htp.structure import StructureEnsemble
from enzy_htp.structure.structure_selection import select_stru
from enzy_htp._interface.amber_interface import AmberParameter
from enzy_htp import (
    interface,
)

ai = interface.amber

# I. Run sander QMMM single point
def qmmm_single_point_sander( # TODO make it a more general API that is independent of Amber or put this inside amber and called but another API (becasue it seems dont need build_md_step other packages)
        stru_esm: StructureEnsemble,
        cluster_job_config: Dict,
    ):
    """Run a QM/MM single point calculate for each frame in {stru_esm}"""
    # 1. build md step
    qm_lot = QMLevelOfTheory(
        basis_set="def2svp",
        method="PBE0",
    )
    qmmm_spe_step = ai.build_md_step(
        name="qmmm_spe",
        minimize=True,
        length=1, # i.e.: single point energy
        cluster_job_config=cluster_job_config,
        core_type="cpu",
        # qmmm settings
        use_qmmm = True,
        qm_region = select_stru(stru_esm.structure_0, "resi 294"),
        qm_region_charge_spin = (-1, 1),
        qm_level_of_theory  = qm_lot,
        qm_engine = "g16",
        qm_ele_cutoff = 12.0,
        qm_ewald = 0,
        qm_adjust_q = 1,
        qm_adaptive_solvent_type = "single_point", # num_of_solvent, fixed_size
        qm_num_adaptive_solvent = 10,
        qm_region_pdb_path="1",
    )

    # 2. do a special iteration of stru_esm that gives paths of inpcrd?
    jobs = []
    output_eggs = []
    for i, stru in enumerate(stru_esm.structures()):
        qmmm_spe_step.work_dir = f"./QMMM/frame_{i}/"
        fs.safe_mkdir(qmmm_spe_step.work_dir)
        inpcrd_path = f"{qmmm_spe_step.work_dir}frame.inpcrd"
        prmtop_path = f"{qmmm_spe_step.work_dir}frame.prmtop"
        ai.convert_stru_to_inpcrd(stru, inpcrd_path) # TODO create temp path
        ai.convert_top_to_prmtop(stru_esm.topology_source_file, prmtop_path) # TODO create temp path
        params = AmberParameter(
            inpcrd_path=inpcrd_path,
            prmtop_path=prmtop_path,
        )

        job, output_egg = qmmm_spe_step.make_job(params, path_rel_to=qmmm_spe_step.work_dir)
        job.sub_dir = qmmm_spe_step.work_dir
        job.sub_script_path = os.path.abspath(job.sub_script_path)
        jobs.append(job)
        output_eggs.append(output_egg)

    # 3. run the md step
    ClusterJob.wait_to_array_end(jobs, period=20, array_size=2)

def _extract_atomic_charge(chg_path: str) -> List:
    """extract atomic charge from a .chg file output by Multiwfn"""
    result = []
    with open(chg_path) as f:
        lines = f.readlines()
    for line in lines:
        atom_data = line.strip().split()
        atom_data = [atom_data[0]] + [float(i) for i in atom_data[1:]]
        result.append(atom_data)

    return result

def calculate_hirshfeld(target_atom_idx_list: List): # TODO make this an API later redesign I/O
    result = []
    # 1. collect result
    work_dirs = [glob.glob(f"QMMM/frame_{i}")[0] for i in range(100)]
    instr_file = "multiwfn_hirshfeld.in"
    settings = [
        "7",
        "1",
        "1",
        "y",
        "0",
        "q",
        "",
    ]
    fs.write_lines(instr_file, settings)
    for work_dir in work_dirs:
        wfn_file = f"{work_dir}/gau_job.fchk"
    # 2. calculate Hirshfeld charge
        log_path = f"{work_dir}/multiwfn_hirshfeld.out"
        chg_path = f"{work_dir}/multiwfn_hirshfeld.chg"

        interface.multiwfn.run_multiwfn(
            wfn_file=wfn_file,
            instr_file=instr_file,
            log_path=log_path,
            cluster_job_config=None,
        )

        fs.safe_mv("gau_job.chg", chg_path)
    # 3. extract charge from log
        chg_list = _extract_atomic_charge(chg_path)
        result.append([chg_list[idx-1][-1] for idx in target_atom_idx_list])
    return result

def main():
    # I. run QMMM
    stru_esm = ai.load_traj("traj.prmtop", "sampled_100f.mdcrd", "traj.pdb")
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "account" : "yang_lab_csb",
            "partition" : "debug",
            'walltime' : '30:00',
        }
    }
    qmmm_single_point_sander(stru_esm, cluster_job_config)
    # II. run analysis
    result = calculate_hirshfeld([9,10])
    save_obj(result, "chg.pickle")
    with open('chg.csv', 'w', encoding='utf-8', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(result)

if __name__ == "__main__":
    main()
