"""The template for EnzyHTP 2.0 shrapnel main script. It creates
children jobs that each handle part of the given tasks.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-11-11
"""
import datetime
from functools import partial
import glob
import itertools
import os
from pathlib import Path
import sys
import numpy as np
import pickle
import uuid
import signal
from typing import Any, Dict, Iterable, List, Tuple, Callable, Union

from enzy_htp.core.exception import ShrapnelChildError
from enzy_htp.mutation import assign_mutant
from enzy_htp.mutation_class import Mutation, get_involved_mutation, generate_from_mutation_flag
from enzy_htp._config.armer_config import ARMerConfig
from enzy_htp.structure import Structure, PDBParser, Atom, StructureConstraint
import enzy_htp.structure.structure_constraint as stru_cons
from enzy_htp.core.clusters.accre_r9 import AccreR9
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
        mutant_inp: Union[Iterable, str] = task["mutants"]
        ligand_chrg_spin_mapper: Dict = task["ncaa_chrgspin"]
        md_constraints: List[Callable[[Structure], StructureConstraint]] = task["constraints"]
        if not isinstance(mutant_inp, str) and len(mutant_inp) > 0 and isinstance(mutant_inp[0][0], Mutation):
            if len(wt_pdb_list) > 1:
                _LOGGER.warning(
                    "mutants are given as list of mutation objects "
                    "but there is more than one wild type! "
                    "This will not work unless their protein part is identical"
                    )
        for wt_pdb in wt_pdb_list:
            stru = PDBParser().get_structure(wt_pdb)
            if isinstance(mutant_inp, str):
                mutants = assign_mutant(stru, mutant_inp, chain_sync_list=chain_sync_list, chain_index_mapper=chain_index_mapper)
            elif len(mutant_inp) > 0 and isinstance(mutant_inp[0][0], Mutation):
                mutants = mutant_inp
            else:
                _LOGGER.error(f"mutants in the task can only be str or List[List[Mutation]]. Got {mutant_inp}")
                raise TypeError
            for mut in mutants:
                child_tasks.append({
                    "uid" : uuid.uuid4(),
                    "stru" : stru, 
                    "mutant" : mut, 
                    "ncaa_chrgspin" : ligand_chrg_spin_mapper, 
                    "constraints" : md_constraints})
    _LOGGER.info(f"Total {len(child_tasks)} tasks will be distributed into {shrapnel_groups} groups.")
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
        _LOGGER.info("Loading existing child jobs from checkpoint.")
        child_jobs: List[ClusterJob] = load_obj(child_job_checkpoint) # NOTE that this will cause source code change in this script cannot effect the content of those existing jobs. (e.g.: partition etc.)
        new_child_jobs = []
        # 1. analyze remaining child_jobs (recycle old one) may benefit from having a mimo
        for job in child_jobs:
            if not job.job_id: # NOTE after the change @2025/7/4, the job id should be always up to date so long as they are submitted through this script.
                job.retrive_job_id() # BUG multiple submissions of the job causes retriving a wrong job id; consider retrieving using a uid in .log?
            if (not job.is_submitted()) or (not job.is_complete()):
                new_child_jobs.append(job)
        child_jobs = new_child_jobs # this will include failing, pending, running jobs
        # 2. update the checkpoint
        save_obj(child_jobs, child_job_checkpoint)

    # signal handling in case hitting the walltime while waiting for child jobs to finish
    _sig_handler_partial = partial(_sig_handler, child_jobs=child_jobs, child_jobs_checkpoint=child_job_checkpoint)
    signal.signal(signal.SIGUSR1, _sig_handler_partial)
    # signal.signal(signal.SIGTERM, _sig_handler_partial)

    jobs_remain = ClusterJob.wait_to_array_end_plus(child_jobs, shrapnel_check_period, shrapnel_child_array_size)
    # update child jobs
    save_obj(child_jobs, child_job_checkpoint)

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
    from enzy_htp.analysis import ele_field_strength_at_along, ddg_fold_of_mutants, rmsd, binding_energy, bond_dipole, ele_stab_energy_of_bond
    from enzy_htp.quantum import single_point
    from enzy_htp.chemical.level_of_theory import QMLevelOfTheory
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

        task_result_data: Dict = result_dict.setdefault(task_uid, dict())
        task_dir = f"task_{i}" # as this will be executed under the grp_dir
    
        # 1. prepare
        prepared_stru = task_result_data.get("prepared_stru", None)
        if prepared_stru is None:
            prepared_stru = remove_hydrogens(wt_stru, polypeptide_only=True)
            protonate_stru(prepared_stru, protonate_ligand=False)
            result_dict[task_uid]["prepared_stru"] = prepared_stru
            save_obj(result_dict, result_path)

        # 2. mutate
        mutant_stru = task_result_data.get("mutant_stru", None)
        if mutant_stru is None:
            mutant_stru = mutate_stru(prepared_stru, mutant, engine="pymol")
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
                    "leaprc.gaff",
                    "leaprc.water.tip3p",
                ],
            )
            task_gpu_job_config = {
                "cluster" : gpu_job_config["cluster"],
                "res_keywords" : gpu_job_config["res_keywords"] | {"node_cores" : gpu_partition}
            }
            task_constraints = []
            for cons in md_constraints:
                task_constraints.append(cons(mutant_stru))
            new_trajs = equi_md_sampling(
                stru = mutant_stru,
                param_method = param_method,
                cluster_job_config = task_gpu_job_config,
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

        # metrics specific
        for metric_key in ("ef", "rmsd", "mmpbsa", "bond_dipole", "dg_ele", "qm_results"):
            task_result_data.setdefault(metric_key, list())

        for j, replica_esm in enumerate(trajs):
            # EF
            if j < len(task_result_data.get("ef", list())):
                replica_ef = task_result_data["ef"][j]
            else:
                replica_ef = []
                atom_1 = mutant_stru.get("C.218.C1")
                atom_2 = mutant_stru.get("C.218.I1")
                ef_region_pattern = "resi 1-216"
                for traj_stru in replica_esm.structures(remove_solvent=True):
                    # 5. get dEF
                    field_strength = ele_field_strength_at_along(
                        traj_stru, atom_1, atom_2, region_pattern=ef_region_pattern)
                    replica_ef.append(field_strength)
                task_result_data["ef"].append(replica_ef)
                save_obj(result_dict, result_path)

            # RMSD
            if j >= len(task_result_data["rmsd"]):
                replica_rmsd = rmsd(
                    replica_esm, 
                    region_pattern="resi 24+25+35+36+61+62+63+64+82+83+84+85+88+110+111+112+113+126+127+128+129+132+169+170+217+218 and (not elem H)",
                )
                task_result_data["rmsd"].append(replica_rmsd)
                save_obj(result_dict, result_path)

            # MMPBSA
            if j >= len(task_result_data["mmpbsa"]):
                replica_mmpbsa = binding_energy(
                    replica_esm,
                    ligand="resi 218",
                    method="mmpbsa_amber",
                    cluster_job_config=cpu_job_config,
                )
                task_result_data["mmpbsa"].append(replica_mmpbsa)
                save_obj(result_dict, result_path)

            # QM & dipole
            if j >= len(task_result_data["qm_results"]):
                qm_results = single_point(
                    stru=replica_esm,
                    engine="gaussian",
                    method=QMLevelOfTheory( basis_set="def2-tzvp", method="pbe0" ),
                    regions=["resi 217+218"],
                    cluster_job_config=cpu_job_config,
                    job_check_period=60,
                    job_array_size=20,
                    work_dir=f"{task_dir}/QM_SPE/rep_{j}",
                )
                task_result_data["qm_results"].append(qm_results)
                save_obj(result_dict, result_path)
            else:
                qm_results = task_result_data["qm_results"][j]

            if j >= len(task_result_data["bond_dipole"]):
                replica_dipole = []
                for ele_stru in qm_results:
                    this_frame_stru = ele_stru.geometry.topology
                    atom_1 = this_frame_stru.get("C.217.C1")
                    atom_2 = this_frame_stru.get("C.217.I1")

                    # bond dipole
                    dipole = bond_dipole(
                        ele_stru, atom_1, atom_2,
                        work_dir=f"{task_dir}/bond_dipole/rep_{j}"
                    )
                    replica_dipole.append(dipole)
                task_result_data["bond_dipole"].append(replica_dipole)
                save_obj(result_dict, result_path)
            else:
                replica_dipole = task_result_data["bond_dipole"][j]

            #dg_ele
            if j >= len(task_result_data["dg_ele"]):
                replica_dg_ele = []
                for dip, ef in zip(replica_dipole, replica_ef):
                    dg_ele = ele_stab_energy_of_bond(dip[0], ef)
                    replica_dg_ele.append(dg_ele)
                task_result_data["dg_ele"].append(replica_dg_ele)
                save_obj(result_dict, result_path)

    # final san check
    for task in tasks:
        task_uid = task["uid"]
        assert task_uid in result_dict
        # assert "ddg_fold" in result_dict[task_uid]
        assert "trajs" in result_dict[task_uid]
        assert len(result_dict[task_uid]["trajs"]) == md_parallel_runs
        assert "ef" in result_dict[task_uid]
        assert len(result_dict[task_uid]["ef"]) == md_parallel_runs
        assert "rmsd" in result_dict[task_uid]
        assert len(result_dict[task_uid]["rmsd"]) == md_parallel_runs
        assert "mmpbsa" in result_dict[task_uid]
        assert len(result_dict[task_uid]["mmpbsa"]) == md_parallel_runs
        assert "bond_dipole" in result_dict[task_uid]
        assert len(result_dict[task_uid]["bond_dipole"]) == md_parallel_runs
        assert "dg_ele" in result_dict[task_uid]
        assert len(result_dict[task_uid]["dg_ele"]) == md_parallel_runs

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

def _dump_state(child_jobs: list, child_job_checkpoint: str):
    """dump important information when the script is terminated."""
    save_obj(child_jobs, child_job_checkpoint)
    _LOGGER.warning(f"[{datetime.datetime.now():%F %T}] Saved {len(child_jobs)} jobs → {child_job_checkpoint}")

def _sig_handler(signum, frame, *, child_jobs, child_jobs_checkpoint):
    """handling signals"""
    _LOGGER.warning(f"Received signal {signum}, dumping state...")
    _dump_state(child_jobs, child_jobs_checkpoint)
    sys.exit(0)  

def resubmit_child_jobs(child_job_list: list, child_job_checkpoint: str, 
                        shrapnel_check_period: int, shrapnel_child_array_size: int, shrapnel_dir: str):
    """Resubmit a specific list of children jobs.
    This function is handy when you have part of the child jobs failed after the main script finishes.
    Always dump resubmitted child jobs to a new file."""
    child_jobs: List[ClusterJob] = load_obj(child_job_checkpoint) # NOTE that this will cause source code change in this script cannot effect the content of those existing jobs. (e.g.: partition etc.)
    child_jobs_mapper = {i.sub_dir.removeprefix(shrapnel_dir) : i for i in child_jobs}

    new_child_jobs = [child_jobs_mapper[group_name] for group_name in child_job_list]
    new_child_job_checkpoint = Path(child_job_checkpoint).with_suffix(".new.pickle")
    new_child_job_checkpoint = fs.get_valid_temp_name(str(new_child_job_checkpoint))
    
    _sig_handler_partial = partial(_sig_handler, child_jobs=new_child_jobs, child_jobs_checkpoint=new_child_job_checkpoint)
    signal.signal(signal.SIGUSR1, _sig_handler_partial)

    jobs_remain = ClusterJob.wait_to_array_end_plus(new_child_jobs, shrapnel_check_period, shrapnel_child_array_size)
    save_obj(new_child_jobs, new_child_job_checkpoint)

    if jobs_remain:
        _LOGGER.error("some children jobs didn't finish normally. they are:")
        _LOGGER.error("\n".join([j.sub_dir for j in jobs_remain]))
        raise ShrapnelChildError

def detect_child_job_progress(
        md_parallel_runs: int,
        child_dir_list: list = None, 
        child_job_checkpoint: str = "shrapnel_child_jobs.pickle",
        child_result_fname: str = "results.pickle",
        child_task_fname: str = "tasks.pickle",
        print_unfinished_only: bool = False,
        return_unfinished: bool = True):
    """Detect the progress of all children jobs. Return the jobs that is not finished.
    Supply child_dir_list or child_job_checkpoint.
    NOTE only detect for MD finish in this template. You can modify this function to customize it."""
    if not child_dir_list:
        if os.path.exists(child_job_checkpoint):
            child_jobs: List[ClusterJob] = load_obj(child_job_checkpoint)
            child_dir_list = [i.sub_dir for i in child_jobs]
        else:
            _LOGGER.error("please supply {child_dir_list} or {child_job_checkpoint}")
            raise ValueError
    unfinished = []
    for child_dir in child_dir_list:
        task_file = f"{child_dir}/{child_task_fname}"
        result_file = f"{child_dir}/{child_result_fname}"
        tasks = load_obj(task_file)
        total_num = len(tasks)
        try:
            results = load_obj(result_file)
            finish_num = 0
            for task in tasks:
                task_finish = _determine_task_md_finish(task, results, md_parallel_runs)
                finish_num += task_finish
        except (pickle.UnpicklingError, EOFError):
            finish_num = -1

        if not print_unfinished_only:
            print(f"{child_dir}: {finish_num}/{total_num}")
        elif finish_num < total_num:
            print(f"{child_dir}: {finish_num}/{total_num}")
            unfinished.append(child_dir)

    if return_unfinished:
        return unfinished

def _determine_task_md_finish(task: dict, results: dict, md_parallel_runs: int,) -> bool:
    task_id = task["uid"]
    task_result = results[task_id]
    trajs = task_result.get("trajs", list())
    if len(trajs) < md_parallel_runs:
        return False
    else:
        return True

def main():
    """shrapnel-like dir creation and submission for"""
    # region: Input
    tasks = [
        {   
        # "wt" : ["wt/athmt.pdb"],
        # "wt" : ["wt/aclmt.pdb"],
        "wt" : ["wt/athmt-eti.pdb"],
        # "mutants" : "a:1[resi 8+11+15+35+127+128+131+160+202:all not self]",
        # "mutants" : load_obj("mutants_hammers_aclmt.pickle"),
        "mutants" : "a:1[resi 132+169+170:all not self]",
        "constraints" : [
            partial(stru_cons.create_distance_constraint,"B.217.S1", "C.218.C1", 3.5),
            partial(stru_cons.create_angle_constraint,"B.217.S1", "C.218.C1", "C.218.I1", 180.0),            
        ],
        "ncaa_chrgspin" : {"SAH" : (0,1), "ETI" : (0,1)},
        },
    ]
    # temp (will deprocate)
    chain_sync_list = []
    chain_index_mapper = {}
    # ARMer settings
    cluster = AccreR9()
    yanglab_acc_res_keywords = {
        "account" : "yang_lab",
        "partition" : "batch",
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
            "partition" : "batch_gpu",
            "walltime" : "3-00:00:00",}}
    # endregion

    run_shrapnel(
        tasks = tasks,
        chain_sync_list=chain_sync_list, chain_index_mapper=chain_index_mapper,
        # MD
        md_length=50.0,
        md_parallel_runs=1,
        # ARMer
        shrapnel_child_job_config = shrapnel_child_job_config,
        shrapnel_cpujob_config = shrapnel_cpujob_config,
        shrapnel_gpujob_config = shrapnel_gpujob_config,
        # shrapnel
        work_dir = "./shrapnel_athmt_3_site",
        shrapnel_child_array_size = 100,
        shrapnel_groups = 100,
        shrapnel_gpu_partition_mapper = {
            (0, 30) : "nvidia_rtx_a6000:1",
            (31, 100) : "nvidia_geforce_rtx_2080_ti:1",
        },
    )

    # resubmit_child_jobs(
    #     child_job_list=[
    #         "group_11",
    #         "group_13",
    #         "group_15",
    #         "group_16",
    #         "group_17",
    #         "group_18",
    #         "group_19",
    #         "group_20",
    #         "group_21",
    #         "group_22",
    #         "group_23",
    #         "group_24",
    #         "group_25",
    #         "group_26",
    #         "group_28",
    #         "group_29",
    #         "group_30",
    #         "group_31",
    #         "group_32",
    #         "group_35",
    #         "group_38",
    #         "group_39",
    #         "group_40",
    #         "group_42",
    #         "group_43",
    #         "group_44",
    #         "group_45",
    #         "group_46",
    #         "group_47",
    #         "group_48",
    #         "group_49",
    #         "group_50",
    #         "group_52",
    #         "group_54",
    #         "group_56",
    #         "group_57",
    #         "group_58",
    #         "group_59",
    #         "group_60",
    #         "group_61",
    #         "group_62",
    #         "group_65",
    #         "group_66",
    #         "group_67",
    #         "group_68",
    #         "group_69",
    #         "group_70",
    #         "group_83",
    #         "group_84",
    #         "group_85",
    #         "group_86",
    #         "group_87",
    #         "group_89",
    #         "group_90",
    #         "group_91",
    #         "group_93",
    #         "group_95",
    #         "group_96",
    #         "group_97",
    #         "group_98",
    #         "group_99",
    #     ],
    #     child_job_checkpoint="shrapnel_child_jobs.pickle",
    #     shrapnel_check_period = 120,
    #     shrapnel_child_array_size = 50,
    #     shrapnel_dir="./shrapnel_aclmt/",
    # )

    unfinished = (detect_child_job_progress(
        md_parallel_runs=1,
        print_unfinished_only=True,
    ))
    # unfinished_idx = set([int(i.removeprefix("./shrapnel_aclmt/group_")) for i in unfinished])
    # running_idx = {
    #     0 ,1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,12,14,27,33,34,36,37,41,51,53,55,63,64,71,72,73,74,75,76,77,78,79,80,81,82,88,92,94
    # }
    # running_md_idx = {
    #     34,92,27,53,36,2,8,0,1,6,5,3,55,7,4,10,9,12,88,76,37,74,71,73,72,77,75,33,78,79,41,80,64,51,63,81,14,82,94
    # }
    # print(unfinished_idx - running_idx)
    # print(unfinished_idx - running_md_idx)
    # print(running_idx - running_md_idx)
    # print(running_md_idx - running_idx)

if __name__ == "__main__":
    main()
