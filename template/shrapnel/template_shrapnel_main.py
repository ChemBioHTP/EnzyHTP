"""The template for EnzyHTP 2.0 shrapnel main script. It creates
children jobs that each handle part of the given tasks.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-11-11
"""
from functools import partial
import glob
import itertools
import os
from pathlib import Path
import numpy as np
import pickle
import uuid
from typing import Any, Dict, List, Tuple, Callable

from enzy_htp.analysis import d_ele_field_upon_mutation_coarse, ddg_fold_of_mutants
from enzy_htp.core.exception import ShrapnelChildError
from enzy_htp.mutation import assign_mutant
from enzy_htp.mutation_class import Mutation, get_involved_mutation, generate_from_mutation_flag
from enzy_htp._config.armer_config import ARMerConfig
from enzy_htp.structure import Structure, PDBParser, Atom, StructureConstraint
import enzy_htp.structure.structure_constraint as stru_cons
from enzy_htp.core.clusters.accre import Accre
from enzy_htp.core.job_manager import ClusterJob, ClusterInterface
from enzy_htp.core import _LOGGER
from enzy_htp.core.general import save_obj, load_obj, save_func_to_main
import enzy_htp.core.math_helper as mh
import enzy_htp.core.file_system as fs

def run_shrapnel(
        tasks: Dict,
        chain_sync_list: List, chain_index_mapper: Dict,
        md_length: float,
        md_parallel_runs: int,
        shrapnel_child_job_config: Dict,
        shrapnel_cpujob_config: Dict,
        shrapnel_gpujob_config: Dict,
        shrapnel_check_period: int = 120,
        shrapnel_child_array_size: int = 50,
        shrapnel_groups: int = 100,
        shrapnel_gpu_partition_mapper: Dict = None,
        ncaa_param_lib_path: str = None,
        out_ref_path_1: str = "result_metrics.pickle",
        child_job_checkpoint: str = "shrapnel_child_jobs.pickle",
        work_dir: str = "./shrapnel",
        # fname under work_dir or child_dir
        child_result_fname = "results.pickle",
        child_task_fname = "tasks.pickle",
        child_main_fname = "child_main.py",
        kwargs_fname = "child_main_kwargs.pickle",
    ) -> List[List[Mutation]]:
    """the main function for running a shrapnel workflow.
    The workflow will:
    1. combine all {wt} and {mutants}
    2. create a child main script
    3. divid all tasks (each is the combination of a "wt" and a "mutant")
       evenly into the {shrapnel_groups}. Each group will contain a pickle
       file with all the information about the task
    
    Args:
        tasks: Dict
            describe all tasks for this shrapnel run. The format is:
            [task, ...]
            "task" in this list:
                {   
                "wt" : [wt_path, ...],
                "mutants" : mutant_pattern,
                "constraints" : [md_constraint, ...],
                "ncaa_chrgspin" : {"RES" : (0,1),...},
                },
            In detail, each key of this dict means:
                wt: List[str]
                    path of wt pdb files, if multiple path exists, they are
                    combined with each mutants.
                mutants: str
                    a EnzyHTP mutation pattern
                ncaa_chrgspin: Dict
                    NCAA chrgspin assignment as a mapper
                constraints: List[Callable[[Structure], StructureConstraint]]
                    constraints in MD as partial functions
        md_length: float
            length of MD
        md_parallel_runs: int
            replica runs of MD
        shrapnel_child_job_config: Dict
            job config for each children jobs of shrapnel
        shrapnel_cpujob_config: Dict
            job config for each cpu jobs of shrapnel
        shrapnel_gpujob_config: Dict
            job config for each gpu jobs of shrapnel
        shrapnel_check_period: int = 120
            job check period for the array of all children jobs
        shrapnel_child_array_size: int = 50
            job array size for the array of all children jobs
        shrapnel_groups: int = 100
            num of groups of shrapnel (i.e., total num of child jobs)
        shrapnel_gpu_partition_mapper: Dict = None
            a mapper that allows assignment of gpu partition names for each children
        ncaa_param_lib_path: str = None
            the path of the ncaa_lib.
            (default: f"{work_dir}/ncaa_lib")
        out_ref_path_1: str = "result_metrics.pickle"
            reference result file
        child_job_checkpoint: str = "shrapnel_child_jobs.pickle"
            check point file that contain all children job objects
        work_dir: str = "./shrapnel"
            the working dir
        child_result_fname = "result.pickle"
            result fname relative to the children dir (i.e., f"{work_dir}/group_{i}/")
        task_fname = "tasks.pickle"
            tasks fname relative to the children dir (i.e., f"{work_dir}/group_{i}/")
        child_main_fname = "child_main.py"
            main script fname relative to the work dir
        kwargs_fname = "child_main_kwargs.pickle"
            fname of the kwargs for the main script relative to the work dir

    Return:
        the mutant space after screening."""
    # 0. init
    if ncaa_param_lib_path is None:
        ncaa_param_lib_path = f"{work_dir}/ncaa_lib"
    # 1. compile child tasks
    child_tasks = []
    for task in tasks:
        wt_pdb_list: str = task["wt"]
        mutant_pattern: str = task["mutants"]
        ligand_chrg_spin_mapper: Dict = task["ncaa_chrgspin"]
        md_constraints: List[Callable[[Structure], StructureConstraint]] = task["constraints"]
        for wt_pdb in wt_pdb_list:
            stru = PDBParser().get_structure(wt_pdb)
            mutants = assign_mutant(stru, mutant_pattern, chain_sync_list=chain_sync_list, chain_index_mapper=chain_index_mapper)
            for mut in mutants:
                child_tasks.append({
                    "uid" : uuid.uuid4(),
                    "stru" : stru, 
                    "mutant" : mut, 
                    "ncaa_chrgspin" : ligand_chrg_spin_mapper, 
                    "constraints" : md_constraints})

    num_task_each_grp = mh.calc_average_task_num(len(child_tasks), shrapnel_groups)
    child_main_path = f"{work_dir}/{child_main_fname}"
    kwargs_file = f"{work_dir}/{kwargs_fname}"

    if not Path(child_job_checkpoint).exists():
        fs.safe_mkdir(work_dir)
        # 2. make child main script
        # make kwargs
        child_main_kwargs = {
            "ncaa_param_lib_path" : os.path.abspath(ncaa_param_lib_path),
            "cpu_job_config" : shrapnel_cpujob_config,
            "gpu_job_config" : shrapnel_gpujob_config,
            "md_length" : md_length,
            "md_parallel_runs" : md_parallel_runs,
        }
        save_obj(child_main_kwargs, kwargs_file)
        # make script
        save_func_to_main(_child_main, os.path.abspath(kwargs_file), child_main_path)
        # 3. distribute tasks into groups & make jobs
        child_jobs = []
        assigned = 0
        for grp_id, num_task in enumerate(num_task_each_grp):
            if num_task == 0:
                continue
            # make dir
            grp_path = f"{work_dir}/group_{grp_id}"
            fs.safe_mkdir(grp_path)
            # dist mutants
            grp_tasks = child_tasks[assigned:assigned+num_task]
            assigned += num_task
            grp_task_pickle = f"{grp_path}/{child_task_fname}"
            save_obj(grp_tasks, grp_task_pickle)
            # make job
            gpu_partition = _find_gpu_partition(grp_id, shrapnel_gpu_partition_mapper)
            sub_script_path = os.path.abspath(f"{grp_path}/submit_main.sh")
            child_jobs.append(
                _make_child_job(
                    grp_id = grp_id,
                    cluster_job_config = shrapnel_child_job_config,
                    child_main_path = os.path.abspath(child_main_path),
                    task_fname = child_task_fname,
                    child_result_fname = child_result_fname,
                    gpu_partition = gpu_partition,
                    sub_dir = grp_path,
                    sub_script_path = sub_script_path,
                )
            )

        # 3. submit jobs
        save_obj(child_jobs, child_job_checkpoint)
    else:
        child_jobs: List[ClusterJob] = load_obj(child_job_checkpoint) # NOTE that this will cause source code change in this script cannot effect the content of those existing jobs. (e.g.: partition etc.)
        new_child_jobs = []
        # 1. analyze remaining child_jobs (recycle old one) may benefit from having a mimo
        for job in child_jobs:
            job.retrive_job_id() # BUG multiple submissions of the job causes retriving a wrong job id; consider retrieving using a uid in .log?
            if (not job.is_submitted()) or (not job.is_complete()):
                new_child_jobs.append(job)
        child_jobs = new_child_jobs # this will include failing, pending, running jobs
        # 2. update the checkpoint
        save_obj(child_jobs, child_job_checkpoint)

    jobs_remain = ClusterJob.wait_to_array_end_plus(child_jobs, shrapnel_check_period, shrapnel_child_array_size)

    if jobs_remain:
        _LOGGER.error("some children jobs didn't finish normally. they are:")
        _LOGGER.error("\n".join([j.sub_dir for j in jobs_remain]))
        raise ShrapnelChildError

    # 4. summarize result
    result_metrics = {}
    for grp_id, num_task in enumerate(num_task_each_grp):
        grp_path = f"{work_dir}/group_{grp_id}"
        result_path = f"{grp_path}/{child_result_fname}"
        child_result = load_obj(result_path)
        result_metrics.update(child_result)
    save_obj(result_metrics, out_ref_path_1)

def _child_main(
        sys_argvs,
        # we dont do type hinting here as they will not be imported at the time of func def in the saved main script
        ncaa_param_lib_path: str,
        cpu_job_config,
        gpu_job_config,
        md_length: float,
        md_parallel_runs: int,
    ):
    """the child main script of the shrapnel treatment
    This will be the content of the main function with
    `-t task_file -p gpu_partition -o result` from cmdline.
    In task_file:
        [
            {"stru" : stru, 
            "mutant" : mut, 
            "ncaa_chrgspin" : ligand_chrg_spin_mapper, 
            "constraints" : md_constraints},
            ...
        ]
    """
    # import section (required by save_func_to_main)
    import os
    from typing import List, Dict, Callable
    from collections import defaultdict

    from enzy_htp.preparation import protonate_stru, remove_hydrogens
    from enzy_htp.mutation import mutate_stru
    from enzy_htp.geometry import equi_md_sampling
    from enzy_htp.analysis import ele_field_strength_at_along, ddg_fold_of_mutants
    from enzy_htp import interface
    from enzy_htp.mutation_class import Mutation
    from enzy_htp.structure import StructureConstraint, Structure, Atom, StructureEnsemble
    from enzy_htp.core.general import load_obj, save_obj

    # cmd inp
    for i, arg in enumerate(sys_argvs):
        if arg == "-t":
            tasks: List[Dict] = load_obj(sys_argvs[i+1])
        if arg == "-p":
            gpu_partition: str = sys_argvs[i+1]
        if arg == "-o":
            result_path: str = sys_argvs[i+1]
    # type hinting
    cpu_job_config: Dict
    gpu_job_config: Dict

    # re-run
    result_dict = defaultdict(dict)
    if os.path.exists(result_path):
        result_dict: Dict = load_obj(result_path)

    # 1. get ddg fold
    # todo_task = []
    # for task in tasks:
    #     task_uid = task["uid"]
    #     mut_result_data: Dict = result_dict.get(task_uid, dict())
    #     if "ddg_fold" not in mut_result_data: # skip in re-run
    #         todo_task.append(task)
    # if todo_task:
    #     grouped_tasks = # make a child task a class and a method that group tasks by same WT
    #     ddg_results = ddg_fold_of_mutants(
    #         wt_stru,
    #         todo_mut,
    #         num_iter = 10,
    #         cluster_job_config = cpu_job_config,
    #         relax_cluster_job_config = cpu_job_config,
    #     )
    #     for k, v in ddg_results.items():
    #         result_dict[k]["ddg_fold"] = v
    #     # save
    #     save_obj(result_dict, result_path)

    for i, task in enumerate(tasks):
        task_uid = task["uid"]
        wt_stru: Structure = task["stru"]
        mutant: List[Mutation] = task["mutant"]
        ligand_chrg_spin_mapper: Dict = task["ncaa_chrgspin"]
        md_constraints: List[Callable[[Structure], StructureConstraint]] = task["constraints"]

        task_result_data: Dict = result_dict.get(task_uid, dict())
        task_dir = f"task_{i}" # as this will be executed under the grp_dir
    
        # 1. prepare
        prepared_stru = task_result_data.get("prepared_stru", None)
        if prepared_stru is None:
            remove_hydrogens(wt_stru, polypeptide_only=True)
            protonate_stru(wt_stru, protonate_ligand=False)
            result_dict[task_uid]["prepared_stru"] = prepared_stru
            save_obj(result_dict, result_path)

        # 2. mutate
        mutant_stru = task_result_data.get("mutant_stru", None)
        if mutant_stru is None:
            mutant_stru = mutate_stru(wt_stru, mutant, engine="pymol")
            mutant_stru.assign_ncaa_chargespin(ligand_chrg_spin_mapper)
            remove_hydrogens(mutant_stru, polypeptide_only=True)
            protonate_stru(mutant_stru, protonate_ligand=False)
            result_dict[task_uid]["mutant_stru"] = mutant_stru
            # save
            save_obj(result_dict, result_path)
        
        # 3. MD
        # TODO check and stop any previous MD jobs (workflow specific)
        trajs: List[StructureEnsemble] = task_result_data.get("trajs", list())
        runs_left = md_parallel_runs - len(trajs)
        if runs_left > 0:
            param_method = interface.amber.build_md_parameterizer(
                ncaa_param_lib_path=ncaa_param_lib_path,
                force_fields=[
                    "leaprc.protein.ff14SB",
                    "leaprc.gaff2",
                    "leaprc.water.tip3p",
                ],
            )
            gpu_job_config = {
                "cluster" : gpu_job_config["cluster"],
                "res_keywords" : gpu_job_config["res_keywords"] | {"partition" : gpu_partition}
            }
            task_constraints = []
            for cons in md_constraints:
                task_constraints.append(cons(mutant_stru))
            new_trajs = equi_md_sampling(
                stru = mutant_stru,
                param_method = param_method,
                cluster_job_config = gpu_job_config,
                prod_constrain=task_constraints,
                prod_time=md_length,
                record_period=md_length*0.01,
                work_dir=f"{task_dir}/MD/",
                parallel_runs=runs_left,
            )
            trajs.extend(new_trajs)
            # save
            result_dict[task_uid]["trajs"] = trajs
            save_obj(result_dict, result_path)

    # final san check
    for task in tasks:
        task_uid = task["uid"]
        assert task_uid in result_dict
        assert "ddg_fold" in result_dict[task_uid]
        assert "trajs" in result_dict[task_uid]

def _make_child_job(
        grp_id: int,
        cluster_job_config: Dict,
        child_main_path: str,
        task_fname: str,
        child_result_fname: str,
        gpu_partition: str,
        sub_dir: str,
        sub_script_path: str,
        ) -> ClusterJob:
    """make the ClusterJob for children runs
    under the shrapnel framework"""
    cluster = cluster_job_config["cluster"]
    res_keywords = cluster_job_config["res_keywords"]

    cmd = f"python -u {child_main_path} -t {task_fname} -p {gpu_partition} -o {child_result_fname} > {child_main_path}.grp{grp_id}.out 2>&1"
    enzyhtp_main_env = cluster.ENZYHTP_MAIN_ENV["CPU"]
    final_res_keywords = ARMerConfig.SINGLE_CPU_RES | {
        'job_name' : f'shrapnel_child_{grp_id}',
        'mem_per_core' : '10G',
        'walltime' : '10-00:00:00',
    } | res_keywords

    job = ClusterJob.config_job(
        commands=cmd,
        cluster=cluster,
        env_settings=enzyhtp_main_env,
        res_keywords=final_res_keywords,
        sub_dir=sub_dir,
        sub_script_path=sub_script_path,
    )
    return job

def _find_gpu_partition(grp_id, partition_mapper) -> str:
    """determine which partition the group should use
    based on the id and the mapper"""
    for (l, h), partition in partition_mapper.items():
        if l <= grp_id <= h:
            return partition

def main():
    """shrapnel-like dir creation and submission for"""
    # region: Input
    tasks = [
        {   
        "wt" : ["wt/test_1.pdb", "wt/test_2.pdb"],
        "mutants" : "WT",
        "constraints" : [
            partial(stru_cons.create_distance_constraint,"C.282.N", "B.281.C47", 2.8),
            partial(stru_cons.create_angle_constraint,"C.282.N", "B.281.C47", "B.281.S7", 180.0),            
            ],
        "ncaa_chrgspin" : {"SAM" : (1,1), "SAI" : (1,1)},
        },
        {   
        "wt" : ["wt/test_3.pdb",],
        "mutants" : "WT, {X###Y}",
        "constraints" : [],
        "ncaa_chrgspin" : {"SAM" : (1,1), "SAI" : (1,1)},
        },
    ]
    # temp (will deprocate)
    chain_sync_list = []
    chain_index_mapper = {}
    # ARMer settings
    cluster = Accre()
    yanglab_acc_res_keywords = {
        "account" : "yang_lab_csb",
        "partition" : "production",
        'walltime' : '10:00:00',
        }
    shrapnel_child_job_config = {
        "cluster" : cluster,
        "res_keywords" : yanglab_acc_res_keywords | {"walltime" : "10-00:00:00"}}
    shrapnel_cpujob_config = {
        "cluster" : cluster,
        "res_keywords" : yanglab_acc_res_keywords | {"walltime" : "2-00:00:00"},}
    shrapnel_gpujob_config = {
        "cluster" : cluster,
        "res_keywords" : yanglab_acc_res_keywords | {
            "account" : "csb_gpu_acc",
            "walltime" : "5-00:00:00",}}
    # endregion

    run_shrapnel(
        tasks = tasks,
        chain_sync_list=chain_sync_list, chain_index_mapper=chain_index_mapper,
        # MD
        md_length=100.0,
        md_parallel_runs=1,
        # ARMer
        shrapnel_child_job_config = shrapnel_child_job_config,
        shrapnel_cpujob_config = shrapnel_cpujob_config,
        shrapnel_gpujob_config = shrapnel_gpujob_config,
        # shrapnel
        shrapnel_child_array_size = 50,
        shrapnel_groups = 100,
        shrapnel_gpu_partition_mapper = {
            (0, 10) : "a6000x4",
            (11, 20) : "turing",
            (21, 30) : "pascal",
            (31, 40) : "a6000x4",
            (41, 50) : "turing",
            (51, 60) : "pascal",
            (61, 70) : "a6000x4",
            (71, 80) : "turing",
            (81, 90) : "pascal",
            (91, 100) : "a6000x4",
        },
    )

if __name__ == "__main__":
    main()
