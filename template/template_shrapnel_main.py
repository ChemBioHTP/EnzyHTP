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

def decode_mutaflag(mutaflag: str) -> tuple:
    old_resi = mutaflag[0]
    new_resi = mutaflag[-1]
    chain_id = mutaflag[1]
    position = mutaflag[2:-1]
    return old_resi, new_resi ,chain_id, position

def parse_mutation(mutations: List[str]):
    result = []
    for mut in mutations:
        old_resi, new_resi, chain_id, position = decode_mutaflag(mut)
        result.append("".join((old_resi, chain_id, position, new_resi)))
    return result

def gather_output(dir_list: List[str]):
    # if not dir_list:
    #     dir_list = [f"{i}/group_{j}" for i in substrate_list for j in range(30)]
    out_dir = "./"
    for wk_dir in dir_list:
        result_path = f"{wk_dir}/Mutation.dat"
        run_cmd(f"cat {result_path} >> {out_dir}ke_result.dat")

def launch_jobs_single_sub(group_idx_list):
    for i in group_idx_list:
        os.chdir(f"group_{i}")
        run_cmd("sbatch sub_enzy_htp.cmd")
        os.chdir("../")

def assign_partition_single_sub(group_idx_list: list, partition: str):
    """assign ACCRE gpu partition for job dirs"""
    for i in group_idx_list:
        script_path = f"./group_{i}/htp_md.py"
        cmd = "sed -i \"s/{\'account\':\'csb_gpu_acc\'}/{\'account\':\'csb_gpu_acc\', \'partition\':\'"+partition+"\'}/\" "+script_path
        # cmd = "sed -i \"s/{\'account\':\'csb_gpu_acc\', \'partition\':\'maxwell\'}/{\'account\':\'csb_gpu_acc\', \'partition\':\'maxwell\', \'walltime\':\'6-00:00:00\'}/\" "+script_path
        # cmd = "sed -i \"s/\'walltime\':\'6-00:00:00\'/\'walltime\':\'5-00:00:00\'/\" "+script_path
        run_cmd(cmd)

def check_not_running_entry(cluster, work_dir, num_group):
    cluster: Accre
    running_jobs = set()
    job_ids = run_cmd("squeue -u shaoq1|grep EF || true").stdout.strip().split("\n")
    if job_ids == [""]:
        job_ids = []
    else:
        job_ids = [x.strip().split()[0] for x in job_ids]
    for job_id in job_ids:
        submit_line = cluster.get_job_info(job_id, "WorkDir:.100")
        if submit_line.startswith(work_dir):
            submit_dir = submit_line.removeprefix(work_dir)
            submit_dir = submit_dir.split("/")[0]
            running_jobs.add(submit_dir)
    if not running_jobs:
        raise Exception("found no running job! Make sure you want to do this!")
    all_dir = set()
    for j in range(num_group):
        all_dir.add(f"group_{j}")
    result = list(all_dir - running_jobs)
    return result

def check_entry_job_ids(cluster, entry_name):
    cluster: Accre
    lyl_job_ids = []
    job_ids = run_cmd("squeue -u shaoq1|grep EHTP_MD|grep R").stdout.strip().split("\n")
    job_ids = [x.strip().split()[0] for x in job_ids]
    for job_id in job_ids:
        submit_line = cluster.get_job_info(job_id, "WorkDir:.100,")
        submit_dir = submit_line.removeprefix("/gpfs52/data/yang_lab/shaoqz/AcPutr_detect/TOT/1st/PuO/MMPBSA/")
        if entry_name in submit_dir:
            lyl_job_ids.append((submit_dir, job_id))

    return lyl_job_ids

def check_dir(dir_list=None, num_group=30):
    result = []
    result_1 = []
    if not dir_list:
        dir_list = [f"group_{j}" for j in range(num_group)]

    for wk_dir in dir_list:
        # if not glob(f"{wk_dir}/mutation*/MD/prod.nc"):
        #     result_1.append(wk_dir)
        # for prod_out_path in glob(f"{wk_dir}/mutation*/MD/prod.out"):
        #     with open(prod_out_path) as f:
        #         if "   5.  TIMINGS" not in f.read():
        #             result_1.append(prod_out_path.removesuffix("/MD/prod.out"))
        if glob(f"{wk_dir}/Mutation.dat"):
            result_1.append(wk_dir)
        # with open(f"{wk_dir}/mmpbsa_gen.py") as f:
        #     if "cycle=" in f.read():
        #         result_1.append(wk_dir)
        # check_files = glob(f"{wk_dir}/slurm*out")
        # q_list = [ 
        #     "FileNotFoundError",
        #     # "Error: an illegal memory access was encountered launching kernel kNLSkinTest", 
        #     # "cudaMemcpy GpuBuffer::Download failed an illegal instruction was encountered",
        #     ]
        # if check_content_of_file_list(check_files, q_list):
        #     result.append(wk_dir)
    # result = set(result_1) - set(result)
    return result_1

def check_content_of_file_list(f_path_list: list, query: list) -> bool:
    detect_flag = 0
    for f_path in f_path_list:
        with open(f_path) as f:
            for q in query:
                if q in f.read():
                    detect_flag = 1
    return detect_flag

def resubmit_job(sub_dirs):
    current_dir = os.getcwd() # in case relative path is using
    for sub_d in sub_dirs:
        os.chdir(sub_d)
        print(f"working on {sub_d}")
        run_cmd("sbatch sub_enzy_htp.cmd")
        os.chdir(current_dir)

def revise_job(sub_dirs):
    for sub_d in sub_dirs:
        script_path = f"{sub_d}/htp_md.py"
        # cmd = f"sed -i 's/PDBMin(/PDBMin(cycle=100000,/' {script_path}"
        # cmd = f"sed -i 's/cycle=100000/cycle=200000/' {script_path}"
        # cmd = f"sed -i 's/period=30/period=120/' {script_path}"
        cmd = f"cp sub_enzy_htp.cmd {sub_d}"
        # cmd = f"sed -i 's/import write_data/import check_complete_metric_run, write_data/' {script_path}" # for rerun
        # cmd_1 = f"""sed -i '20a \ \ \ \ \ \ \ \ if check_complete_metric_run(mut, data_output_path):\\
        #     continue' {script_path}"""
        # script_path = f"{sub_d}/sub_enzy_htp.cmd"
        # cmd = f"sed -i 's/3-00/6-00/' {script_path}"
        print(f"working on {sub_d}")
        run_cmd(cmd)
        # run_cmd(cmd_1)

def find_md_job_mutation_substrate(job_dir, submit_line):
    substrate = job_dir.split("/")[0]
    muta_tag = re.findall(r'.+\/mutation_(.*)\/MD\/.*', submit_line)[0]
    muta_tags = repr(list(muta_tag.split("_")))
    return substrate, muta_tags

def find_md_job_mutation_single_substrate(submit_line):
    muta_tag = re.findall(r'.+\/mutation_(.*)\/MD\/.*', submit_line)[0]
    muta_tags = repr(list(muta_tag.split("_")))
    sub_dir = submit_line.split("/")[1]
    return sub_dir, muta_tags

def find_qm_job_mutation_single_substrate(submit_line):
    muta_tag = re.findall(r'.+\/mutation_(.*)\/QM_cluster\/.*', submit_line)[0]
    muta_tags = repr(list(muta_tag.split("_")))
    sub_dir = submit_line.split("/")[1]
    return sub_dir, muta_tags

def check_complete_md_job(md_dir: str):
    """check if the md job in the md_dir is complete"""
    prod_out = f"{md_dir}prod.out"
    if os.path.exists(prod_out):
        with open(prod_out) as f:
            if "   5.  TIMINGS" in f:
                return True
    else:
        return False

def cat_md_run(cluster, not_running_list):
    cat_script = "../cat_qm.py"
    sub_script = "../sub_cat_qm.sh"
    for job_dir in not_running_list:
        with open(f"{job_dir}/submitted_job_ids.log") as f:
            last_job = f.readlines()[-1].strip()
            # print(f"{job_dir} last job: {last_job}")
            job_id, submitline = last_job.split()
            job_state = cluster.get_job_info(job_id, "State")
            print(f"{job_dir} last state: {job_state}")
        if job_state == "COMPLETED":
            if "PDBMD_GPU" in last_job:
                sub_dir, muta_tags = find_md_job_mutation_single_substrate(submitline)
                sub_dir = f"{job_dir}/{sub_dir}"
                if check_complete_md_job(f"{sub_dir}/MD/"):
                    print(f'{job_dir} last MD job finished: cating QM...')
                    # work on cat the job
                    sub_script_path = shutil.copy(sub_script, sub_dir)
                    script_path = shutil.copy(cat_script, sub_dir)
                    run_cmd(f'sed -i "s/XXX/{muta_tags}/" {script_path}')
                    run_cmd(f'cd {sub_dir} && sbatch .{sub_script_path.removeprefix(sub_dir)} && cd ../../')
                else:
                    print(f'{job_dir} last MD job appears to finish but not: resubmitting...')
                    resubmit_job([job_dir])

            if "PDBMin" in last_job:
                print(f'{job_dir} last Min job finished: resubmitting...')
                resubmit_job([job_dir])
            
            if "qm_cluster" in last_job:
                print(f'{job_dir} breakpoint during QM: cating QM...')
                # work on cat the job
                sub_dir, muta_tags = find_qm_job_mutation_single_substrate(submitline)
                sub_dir = f"{job_dir}/{sub_dir}"
                sub_script_path = shutil.copy(sub_script, sub_dir)
                script_path = shutil.copy(cat_script, sub_dir)
                run_cmd(f'sed -i "s/XXX/{muta_tags}/" {script_path}')
                run_cmd(f'cd {sub_dir} && sbatch .{sub_script_path.removeprefix(sub_dir)} && cd ../../')

            # if "MMPBSA" in last_job:
            #     # work on cat the job
            #     with open(f"{job_dir}/Mutation.dat") as f:
            #         if len(re.findall("TAG", f.read())) == len(glob(f"{job_dir}/mutation*/")):
            #             print(f"{job_dir} all current MD entry have finished data. Resubmitting")
            #             run_cmd(f'cd {job_dir} && sbatch sub_enzy_htp.cmd && cd ../../')
            #         else:
            #             print(f"--> work here {job_dir}: job dies before MMPBSA finishes. cp result in.")
            pass
        if job_state == "PENDING":
            if "MMPBSA" in last_job:
                pass
            else:
                print(f"resubmitting {job_dir}. killing pending md job.")
                cluster.kill_job(job_id)
                resubmit_job([job_dir])
        if job_state == "FAILED":
            failed_reason = run_cmd(f"cat {job_dir}/slurm-{job_id}.out|grep -E 'Error|error|ERROR' || true").stdout
            node_list = cluster.get_job_info(job_id, "NodeList")
            print(f"{job_dir} last job failed on {node_list}: {failed_reason}")
        
        if job_state in ["TIMEOUT", "CANCELLED"]:
            if "qm_cluster" in last_job:
                print("<---- here")
            print(f'{job_dir} last job {job_state}: resubmitting...')
            resubmit_job([job_dir])

def not_completed(job_dirs):
    result = []
    for job_dir in job_dirs:
        total_num = None
        current_num = 0
        data_path = f"{job_dir}/Mutation.dat"
        script_path = f"{job_dir}/htp_md.py"
        with open(script_path) as f:
            mutation_list = re.search("mutants = (.*)", f.read()).group(1)
        total_num = len(eval(mutation_list))
        if os.path.exists(data_path):
            with open(data_path) as f:
                current_num = len(re.findall("TAG", f.read()))
            if current_num == total_num:
                continue
        result.append(job_dir)#, total_num, current_num))
    return result

def do_completed(job_dirs):
    result = []
    for job_dir in job_dirs:
        total_num = None
        current_num = 0
        data_path = f"{job_dir}/Mutation.dat"
        script_path = f"{job_dir}/htp_md.py"
        with open(script_path) as f:
            mutation_list = re.search("mutants = (.*)", f.read()).group(1)
        total_num = len(eval(mutation_list))
        if os.path.exists(data_path):
            with open(data_path) as f:
                current_num = len(re.findall("TAG", f.read()))
            if current_num == total_num:
                result.append(job_dir)
    return result

def generate_input_dirs_single_sub(starting_pdb: str, mutations: List[List[str]], num_group: int= 30, work_dir: str=None):
    """generate working dirs for a list of mutations"""
    if work_dir is None:
        work_dir = f"./{starting_pdb.removesuffix('.pdb')}/"
    mkdir(work_dir)

    num_mut_each_group: list = assign_task_num(len(mutations), num_group)
    already_assigned = 0
    for i, mut_num in enumerate(num_mut_each_group):
        dir_name = f"{work_dir}group_{i}"
        mkdir(dir_name)
        # deal with ending group
        group_mutants = mutations[already_assigned:already_assigned+mut_num]
        already_assigned += mut_num
        group_mutants = [parse_mutation(list(x), resi_idx_mapper['ke']) if len(x) != 0 else [] for x in group_mutants]
        # move source pdb
        shutil.copy(starting_pdb, dir_name)
        # move script
        script_path = shutil.copy("./htp_md.py", dir_name)
        run(f'sed -i "s/XXX/{repr(group_mutants)}/" {script_path}', check=True, text=True, shell=True, capture_output=True)
        run(f'sed -i "s/YYY/{starting_pdb}/" {script_path}', check=True, text=True, shell=True, capture_output=True)
        # move submission script
        shutil.copy("./sub_enzy_htp.cmd", dir_name)
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
    num_group = 50
    # == dir gen ==
    with open("mutant_list_2nd.pickle", "rb") as f:
        mutants_list = pickle.load(f)
    for idx, mutants in enumerate(mutants_list):
        mutants = [i for i,j,k in mutants]
        generate_input_dirs_single_sub("KE_07_R7_2_S.pdb", mutants, num_group, f"./sample_{idx}/")

    # == revise and submission ==
    # revise_job([f"group_{j}" for j in range(50)])
    # revise_job([f"sample_{i}/group_{j}" for i in range(5) for j in range(50)])
    # assign_partition_single_sub(range(40,50), "turing")
    # launch_jobs_single_sub(range(0,50))

    # === job maintenance ===
    # cluster= Accre()
    # not_running_job = not_completed(check_not_running_entry(cluster, f"{os.path.abspath(os.curdir)}/", 50))
    # # resubmit_job(not_running_job)
    # # cat_md_run(cluster, not_running_job)#["spe/group_19"])
    # # resubmit_job(not_completed([f"group_{i}" for i in range(50)]))
    # # print(*not_completed([f"sample_{j}/group_{i}" for i in range(50) for j in range(2)]), sep="\n")
    # print(*not_running_job, sep="\n")
    
    # === collect results ===
    # gather_output(check_dir(num_group=40))

if __name__ == "__main__":
    main()
