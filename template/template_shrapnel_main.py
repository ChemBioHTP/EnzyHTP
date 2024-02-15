'''The template for the old EnzyHTP shrapnel main script:
This script creates multiple sub-directories that each contains an HPC-job submission script
corresponding to a workflow that calculate properties for mutants based on QM and MM.
In this way, for example 100 mutants can be split into groups of 5. Up to 20 HPC-jobs can be
submit simultaneously to maximize the efficiency.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-8-8'''
import re
from typing import List
import os
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

def assign_partition(wt_list: List[str], group_idx_list: List[int], partition: str, script_rel_path: str):
    """assign ACCRE gpu partition for job dirs"""
    for wt in wt_list:
        wt_dir = wt.removesuffix(".pdb")
        for i in group_idx_list:
            script_path = f"{wt_dir}/group_{i}/{script_rel_path}"
            cmd = "sed -i \"s/{\'account\':\'csb_gpu_acc\'}/{\'account\':\'csb_gpu_acc\', \'partition\':\'"+partition+"\'}/\" "+script_path
            run_cmd(cmd)
            cmd = "sed -i -E \"s/\\{\'account\':\'csb_gpu_acc\', *\'partition\':\'.+\'\\}/{\'account\':\'csb_gpu_acc\', \'partition\':\'"+partition+"\'}/\" "+script_path
            run_cmd(cmd)

def submit_jobs(wt_list: List[str], group_idx_list: List[int], sub_script_rel_path: str):
    for wt in wt_list:
        wt_dir = wt.removesuffix(".pdb")
        for i in group_idx_list:
            os.chdir(f"{wt_dir}/group_{i}")
            run_cmd(f"sbatch {sub_script_rel_path}")
            os.chdir("../../")

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

def generate_sub_wkdirs(wt_pdb_list: List[str],
                        mutants: List[List[str]],
                        child_script: str,
                        submission_script: str,
                        num_group: int = 30,
                        work_dir: str = None):
    """Generate sub- work dirs for a list of mutations and WTs. The {mutants} will be split
    based on {num_group} into each sub- work dirs.
    Args:
        wt_pdb: the path of the WT PDB file
        mutants: a list of mutants each described by a list of mutations described by
            a string in the format of {orig_resi_one_letter_name}{chain_id}{orig_resi_index}{target_resi_one_letter_name}.
            e.g.: SA48N
        child_script: the template child script serve as main script of each group
        submission_script: the submission script of each group
        num_group: number of groups for dividing the mutants"""
    if work_dir is None:
        work_dir = "./shrapnel_htp"
    mkdir(work_dir)
    
    for wt_pdb in wt_pdb_list:
        # prepare workdir
        sub_work_dir = f"{work_dir}/{wt_pdb.removesuffix('.pdb')}/"
        mkdir(sub_work_dir)

        # make groups
        num_mut_each_group: list = assign_task_num(len(mutants), num_group)
        already_assigned = 0
        for i, mut_num in enumerate(num_mut_each_group):
            dir_name = f"{sub_work_dir}group_{i}"
            mkdir(dir_name)

            # assign mutants
            group_mutants = mutants[already_assigned:already_assigned+mut_num]
            group_mutants_pickle_path = f"{dir_name}/group_mutants.pickle"
            with open(group_mutants_pickle_path, "wb") as of:
                pickle.dump(group_mutants, of)
            already_assigned += mut_num
            # copy source pdb
            shutil.copy(wt_pdb, dir_name)
            # copy/modify script
            script_path = shutil.copy(child_script, dir_name)
            run_cmd(f'sed -i "s/YYY/{wt_pdb}/" {script_path}')
            # move submission script
            group_sub_script_path = shutil.copy(submission_script, dir_name)
            run_cmd(f'sed -i "s/XXX/{i}_{wt_pdb.removesuffix(".pdb")}/" {group_sub_script_path}')
            # move ./ligands and ./maas
            mkdir(f"{dir_name}/ligands")
            run_cmd(f'cp -r ./ligands/* {dir_name}/ligands')
            mkdir(f"{dir_name}/maas")
            run_cmd(f'cp -r ./maas/* {dir_name}/maas')

def assign_task_num(num_of_task: int, num_of_worker: int) -> List[int]: # TODO probably can also use itertools.cycle
    """assign task number for each worker based on {num_of_task} and {num_of_worker}"""
    task_each_worker = num_of_task // num_of_worker
    remaining_tasks = num_of_task % num_of_worker
    result = [task_each_worker for i in range(num_of_worker)]
    for i in range(remaining_tasks):
        result[i] += 1
    return result

def find_not_running_entry( # Specific to ACCRE
        cluster: Accre,
        job_list_cmd: str,
        target_dir: str,
        wt_list: List[str],
        num_group: int,
        if_print_sline: bool=False):
    """find all the not running EnzyHTP shrapnel sub-workdir
    Args:
        cluster: the object that injects information about the current HPC.
        job_list_cmd: the command that find a set of initial jobs that contains all jobs of the target shrapnel project.
                    properly setting this will reduce the time cost of this function
                    (example: squeue -u username)
        target_dir: the absolute path of the target parent directory of all the sub-workdirs.
                    (used to distinguish the target project from other shrapnel projects.)
        wt_list: the list of wt pdb paths that are used to determine the structure of sub-workdirs
        num_group: the number of subdir groups under each wt
        if_print_sline: a debug option that turns on/off the print of the job info line"""
    running_jobs = set()
    # 1. find all running job ids
    job_ids = run_cmd(f"{job_list_cmd} || true").stdout.strip().split("\n")
    if job_ids == [""]:
        job_ids = []
    else:
        job_ids = [x.strip().split()[0] for x in job_ids]

    # 2. determine the working dir for each running job
    for job_id in job_ids:
        submit_line = cluster.get_job_info(job_id, "WorkDir:.100")
        if if_print_sline:
            print(f"{job_id} {submit_line}")
        if submit_line.startswith(target_dir):
            submit_dir = submit_line.removeprefix(target_dir)
            running_jobs.add(submit_dir)

    if not running_jobs:
        raise Exception("found no running job! Make sure you want to do this!")
    
    # 3. deduce not running jobs
    all_dir = set()
    for wt in wt_list:
        wt_dir = wt.removesuffix(".pdb")
        for idx in range(num_group):
            all_dir.add(f"{wt_dir}/group_{idx}")
    result = list(all_dir - running_jobs)

    print(f"Total: {len(result)}")

    return result

def resubmit_job(sub_dirs: str, submission_script):
    """submit job based on a job dir. Mainly used for resubmission"""
    current_dir = os.getcwd() # in case relative path is using
    for sub_d in sub_dirs:
        os.chdir(sub_d)
        print(f"working on {sub_d}")
        run_cmd(f"sbatch {submission_script}")
        os.chdir(current_dir)

def main(): # TODO update this to master
    num_group = 20 # the number of groups
    # child_script="template_child_w_rcons_main.py"
    child_script="template_child_main.py"
    submission_script="template_hpc_submission.sh"
    data_rel_path="mutant_property.dat"

    # == generate sub-directories ==
    wt_list = ["3cfr-rlp-pea_ah.pdb", "3cfr-slp-pea_ah.pdb"]

    # with open("mutants_ata.pickle", "rb") as f:
    #     mutants = pickle.load(f)
    # generate_sub_wkdirs(wt_list,
    #                     mutants,
    #                     child_script,
    #                     submission_script,
    #                     num_group,
    #                     "./wo_react_cons")

    # == revision and submission ==
    # assign_partition(wt_list, range(0,10), "pascal",
    #                  script_rel_path=child_script)
    # submit_jobs(wt_list, range(5,20),
    #             sub_script_rel_path=submission_script)

    # === maintenance ===
    cluster= Accre()

    # find not running
    not_running_jobs = find_not_running_entry(
        cluster=cluster,
        job_list_cmd="squeue -u shaoq1 | grep 3cfr",
        target_dir=f"{os.path.abspath(os.curdir)}/",
        wt_list=wt_list,
        num_group=num_group,
        if_print_sline=False)
    print(*not_running_jobs, sep="\n")
    # resubmit_job(not_running_jobs, submission_script)

    # # determine complete
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
