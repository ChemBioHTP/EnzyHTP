'''The template for the old EnzyHTP shrapnel main script:
This script creates multiple sub-directories that each contains an HPC-job submission script
corresponding to a workflow that calculate properties for mutants based on QM and MM.
In this way, for example 100 mutants can be split into groups of 5. Up to 20 HPC-jobs can be
submit simultaneously to maximize the efficiency.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-8-8'''
from glob import glob
import re
from typing import Dict, List
import os
from subprocess import run, CalledProcessError
import shutil
import pickle

from helper import mkdir, run_cmd
from core.clusters.accre import Accre

def check_group_w_data(group_list: str, data_rel_path: str):
    "return groups that contain data"
    result = []
    for wk_dir in group_list:
        if os.path.exists(f"{wk_dir}/{data_rel_path}"):
            result.append(wk_dir)
    return result

def gather_output(dir_list: List[str], data_rel_path: str, out_path: str):
    for wk_dir in dir_list:
        result_path = f"{wk_dir}/{data_rel_path}"
        run_cmd(f"cat {result_path} >> {out_path}")

def assign_partition(group_idx_list: List[int], partition: str, script_rel_path: str):
    """assign ACCRE gpu partition for job dirs"""
    for i in group_idx_list:
        script_path = f"./group_{i}/{script_rel_path}"
        cmd = "sed -i \"s/{\'account\':\'csb_gpu_acc\'}/{\'account\':\'csb_gpu_acc\', \'partition\':\'"+partition+"\'}/\" "+script_path
        run_cmd(cmd)

def submit_jobs(group_idx_list: List[int], sub_script_rel_path: str):
    for i in group_idx_list:
        os.chdir(f"group_{i}")
        run_cmd(f"sbatch {sub_script_rel_path}")
        os.chdir("../")

def check_complete_group(group_list: List[str],
                         script_rel_path: str,
                         data_rel_path: str,
                         return_not_complete: bool= False) -> List[str]:
    """return group paths from {group_list} that are completed (or not)"""
    result = []
    for job_dir in group_list:
        current_num = 0
        data_path = f"{job_dir}/{data_rel_path}"
        script_path = f"{job_dir}/{script_rel_path}"
        with open(script_path) as f:
            mutation_list = re.search("mutants = (.*)", f.read()).group(1)
        total_num = len(eval(mutation_list))
        if os.path.exists(data_path):
            with open(data_path) as f:
                current_num = len(re.findall("TAG", f.read()))
            if return_not_complete:
                if current_num != total_num:
                    result.append(job_dir)
            else:
                if current_num == total_num:
                    result.append(job_dir)

    return result

def generate_sub_wkdirs(wt_pdb: str,
                        mutants: List[List[str]],
                        child_script: str,
                        submission_script: str,
                        num_group: int = 30,
                        work_dir: str = None,):
    """Generate sub- work dirs for a list of mutations. The {mutants} will be split
    based on {num_group} into each sub- work dirs.
    Args:
        wt_pdb: the path of the WT PDB file
        mutants: a list of mutants each described by a list of mutations described by
            a string in the format of {orig_resi_one_letter_name}{chain_id}{orig_resi_index}{target_resi_one_letter_name}.
            e.g.: SA48N
        child_script: the template child script serve as main script of each group
        submission_script: the submission script of each group
        num_group: number of groups for dividing the mutants
        work_dir: the workdir that contains the generated sub- workdirs"""
    # prepare workdir
    if work_dir is None:
        work_dir = f"./{wt_pdb.removesuffix('.pdb')}/"
    mkdir(work_dir)

    # make groups
    num_mut_each_group: list = assign_task_num(len(mutants), num_group)
    already_assigned = 0
    for i, mut_num in enumerate(num_mut_each_group):
        dir_name = f"{work_dir}group_{i}"
        mkdir(dir_name)

        # assign mutants
        group_mutants = mutants[already_assigned:already_assigned+mut_num]
        already_assigned += mut_num
        # copy source pdb
        shutil.copy(wt_pdb, dir_name)
        # copy/modify script
        script_path = shutil.copy(child_script, dir_name)
        run_cmd(f'sed -i "s/XXX/{repr(group_mutants)}/" {script_path}')
        run_cmd(f'sed -i "s/YYY/{wt_pdb}/" {script_path}')
        # move submission script
        shutil.copy(submission_script, dir_name)
        # move ./ligands
        mkdir(f"{dir_name}/ligands")
        run_cmd(f'cp -r ./ligands/* {dir_name}/ligands')

def assign_task_num(num_of_task: int, num_of_worker: int) -> List[int]:
    """assign task number for each worker based on {num_of_task} and {num_of_worker}"""
    task_each_worker = num_of_task // num_of_worker
    remaining_tasks = num_of_task % num_of_worker
    result = [task_each_worker for i in range(num_of_worker)]
    for i in range(remaining_tasks):
        result[i] += 1
    return result

def main():
    num_group = 5 # the number of groups
    child_script="template_child_main.py"
    submission_script="template_hpc_submission.sh"
    data_rel_path="Mutation.dat"

    # == generate sub-directories ==
    with open("mutant_list.pickle", "rb") as f:
        mutants = pickle.load(f)
    generate_sub_wkdirs("KE_07_R7_2_S.pdb",
                        mutants,
                        child_script,
                        submission_script,
                        num_group)

    # == revision and submission ==
    # assign_partition(range(3,5), "turing",
    #                  script_rel_path=child_script)
    submit_jobs(range(0,5),
                sub_script_rel_path=submission_script)

    # === maintenance ===
    # cluster= Accre()
    # complete_groups = check_complete_group([f"group_{i}" for i in range(5)],
    #                                        child_script,
    #                                        data_rel_path,
    #                                        return_not_complete=False)
    
    # === collect results ===
    # gather_output(check_group_w_data([f"group_{i}" for i in range(5)], data_rel_path),
    #               data_rel_path,
    #               "./result.dat")

if __name__ == "__main__":
    main()
