"""the main script for the MultiMutDes-EF workflow.
Specialized for the optimization of PuO.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-11
"""
from functools import partial
import glob
import itertools
import os
from pathlib import Path
import numpy as np
import pickle
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

sp = PDBParser()

def ef_optimzer(stru: Structure, atom_1: str, atom_2: str, cutoff: float,
                chain_sync_list: List[tuple],
                chain_index_mapper: Dict[str, int],
                limit_size: bool = True
                ) -> List[List[Mutation]]:
    """optimze the mutant space for maximize the
    change of the internal EF in terms of it's stablization
    of a breaking bond.
    Args:
        stru:
            the structure of the target enzyme-substrate complex
        atom_1, atom_2:
            the target bond defined by 2 atoms (atom_1(-) -> atom_2(+))
            *NOTE has to be (-) -> (+)
        cutoff:
            the cutoff value of the angle.
        chain_sync_list, chain_index_mapper:
            arguments to sync mutants to other non-monomer chains
        limit_size:
            whether limit the size change of mutation
    Returns:
        the mutant space
        """
    if cutoff <= 0:
        _LOGGER.error("cutoff must be a positive value")
        raise ValueError
    cutoff_1 = 180.0 - cutoff
    cutoff_2 = 0.0 + cutoff
    if limit_size:
        size_pattern = " and similar_size_20"
    else:
        size_pattern = ""
    ef_pattern = (f"a:[$ef_hotspot('{atom_1}','{atom_2}', ({cutoff_1},180.0)):charge+{size_pattern},"
                  f" $ef_hotspot('{atom_1}','{atom_2}', (0.0,{cutoff_2})):charge-{size_pattern}]")
    mutants = assign_mutant(stru, ef_pattern,
                            chain_sync_list,
                            chain_index_mapper,)

    return mutants

def coarse_filter(
        stru: Structure,
        mutants: List[List[Mutation]],
        atom_1: str, atom_2: str,
        cluster_job_config: Dict,
        ddg_sample_ranges: List[tuple] = [(float("-inf"),0), (0,10), (10,20), (20,30), (30,40)],
        sample_per_group: int = 300,
        out_ref_path_1: str = "mutant_space_1_coarse_def_ddg_fold.pickle",
        out_ref_path_2: str = "mutant_space_2_coarse_def_ddg_fold.pickle",
        ddg_fold_file: str = None,
    ) -> List[List[List[Mutation]]]:
    """the coarse filter applied on the 1st mutant space ({mutants})
    and assess dEF and ddG_fold. Gives the 2nd mutant space.
    if `ddg_fold_file` is given, ddg_fold values for single mutations will be taken from the file.

    Returns:
        - the 2nd mutant space as a list of sample groups of mutants
        - also save the {mutant : (dEF, ddG_fold)} to a .pickle
          file as a reference.
    Details:
        screening through the 1st mutant space with coarse scores and 
        generate the 2nd mutant space:
        1. stability cut-off: by default 5 groups of 40-30 REU, 30-20 REU, 20-10 REU, 10-0 REU, <0 REU
        2. electrostatic ranking: top ranking {sample_per_group} mutants in each group
    """
    result = []
    result_dicts = []
    # generate {mutant : (dEF, ddG_fold)} dict
    coarse_metrics_dict = _generate_coarse_metrics(stru, mutants, atom_1, atom_2, cluster_job_config, ddg_fold_file)
    save_obj(coarse_metrics_dict, out_ref_path_1)

    # filter mutants
    ddg_groups = [[] for _ in range(len(ddg_sample_ranges))]
    # filter #1
    # stability cut-off: groups of {ddg_sample_ranges}
    for i, (low, high) in enumerate(ddg_sample_ranges):
        ddg_groups[i] = dict(filter(lambda x: low < x[1][1] <= high, coarse_metrics_dict.items()))
    # filter #2
    # electrostatic ranking: top ranking {sample_per_group} mutants in each group
    for ddg_group in ddg_groups:
        sort_by_ef = sorted(ddg_group.items(), key=lambda x: x[1][0], reverse=True)
        result_dicts.append(dict(sort_by_ef[:sample_per_group]))
        result.append([mut for mut, metrics in sort_by_ef])

    save_obj(result_dicts, out_ref_path_2)
    return result

def _generate_coarse_metrics(
        stru: Structure,
        mutants: List[List[Mutation]],
        atom_1: str, atom_2: str,
        cluster_job_config: Dict,
        ddg_fold_file: str,
    ) -> Dict[Tuple[Mutation], Tuple[float, float]]:
    """calculate the dEF and ddG_fold for each mutant
    in mutants"""
    # init
    atom_1: Atom = stru.get(atom_1)
    atom_2: Atom = stru.get(atom_2)
    dp_start = np.array(atom_1.coord)
    dp_end = np.array(atom_2.coord)
    point_of_ef_measure = (dp_start + dp_end)/2
    vec_of_ef_measure = dp_end - dp_start
    result = {}
    for mutant in mutants:
        result[tuple(mutant)] = []

    # calculate dEF
    mutation_ef_mapper = {}
    for mutant in mutants:
        mutant_ef_strength = 0.0 # this also makes WT in [] form be 0.0 
        for mut in mutant:
            if mut not in mutation_ef_mapper:
                mutation_ef_mapper[mut] = d_ele_field_upon_mutation_coarse(
                    stru, mut, point_of_ef_measure, vec_of_ef_measure,
                    unit="MV/cm",
                )
            mutant_ef_strength += mutation_ef_mapper[mut]
        result[tuple(mutant)].append(mutant_ef_strength)
    
    # calculate ddg_fold
    single_point_stability_mapper = _get_single_mutation_stability(
        stru=stru,
        mutant_space=mutants,
        cluster_job_config=cluster_job_config,
        ddg_fold_file=ddg_fold_file,
    )
    for mutant in mutants:
        mutant_stability = 0.0 # this also makes WT in [] form be 0.0 
        for mut in mutant:
            if mut not in single_point_stability_mapper:
                _LOGGER.error(f"{mut} not in single_point mapper!")
                raise Exception
            mutant_stability += single_point_stability_mapper[mut]
        result[tuple(mutant)].append(mutant_stability)

    return result

def _get_single_mutation_stability(
        stru: Structure,
        mutant_space: List[List[Mutation]],
        cluster_job_config: Dict,
        job_check_period: int= 210, # s
        job_array_size: int= 100,
        ddg_fold_file: str= None,
    ) -> Dict[Mutation, float]:
    """calculate ddg_fold for each single point mutations
    involved in the {mutant_space} of {stru}"""
    single_mutations = [[mut] for mut in get_involved_mutation(mutant_space)]
    if ddg_fold_file:
        # 1. extract
        result = load_obj(ddg_fold_file)
        # 2. safety check
        for mut in single_mutations:
            mut = mut[0]
            if mut not in result:
                _LOGGER.error(f"ddg_fold_file missing mutation: {mut} ({ddg_fold_file})")
                raise ValueError
    else:
        # TODO after finish the engine
        ddg_fold_mapper = ddg_fold_of_mutants(
            stru,
            mutant_space=single_mutations,
            method="rosetta_cartesian_ddg",
            cluster_job_config=cluster_job_config,
            job_check_period=job_check_period,
            job_array_size=job_array_size,
            )
        # TODO after finish the engine
        result = {k[0] : v for k, v in ddg_fold_mapper.items()}
    return result

def fine_filter( # TODO summarize logics from here to a general shrapnel function & core checkpoit logics
        stru: Structure,
        mutants: List[List[Mutation]],
        atom_1: str, atom_2: str,
        ligand_chrg_spin_mapper: Dict,
        md_constraints: List[Callable[[Structure], StructureConstraint]],
        md_length: float,
        md_parallel_runs: int,
        ef_region_pattern: str,
        shrapnel_child_job_config: Dict,
        shrapnel_cpujob_config: Dict,
        shrapnel_gpujob_config: Dict,
        shrapnel_check_period: int = 120,
        shrapnel_child_array_size: int = 50,
        shrapnel_groups: int = 100,
        shrapnel_gpu_partition_mapper: Dict = None,
        ncaa_param_lib_path: str = None,
        out_ref_path_1: str = "mutant_space_2_fine_metrics.pickle",
        out_ref_path_2: str = "mutant_space_3_fine_metrics.pickle",
        checkpoint_1: str = "fine_filer_child_jobs.pickle",
        work_dir: str = "./fine_filter",
        # fname under work_dir or child_dir
        child_result_fname = "result.pickle",
        mutant_fname = "mutants.pickle",
        child_main_fname = "child_main.py",
        kwargs_fname = "child_main_kwargs.pickle",
    ) -> List[List[Mutation]]:
    """the fine filter applied on the 2nd mutant space
    and assess 
        dEF(atom_1, atom_2),
        ddG_fold, 
        MMPBSA(ligand), 
        and SPI (ligand, pocket).
    Gives the final mutant space.

    Args:
        stru: Structure
            target wt strutcure
        mutants: List[List[Mutation]]
            target mutant space
        atom_1: str, atom_2: str
            target bond for dEF
        ligand_chrg_spin_mapper: Dict
            NCAA chrgspin assignment
        md_constraints: List[Callable[[Structure], StructureConstraint]]
            constraints in MD
        md_length: float
            length of MD
        md_parallel_runs: int
            replica runs of MD
        ef_region_pattern: str
            the region of the field source charges
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
        out_ref_path_1: str = "mutant_space_2_fine_metrics.pickle"
            reference result file
        out_ref_path_2: str = "mutant_space_3_fine_metrics.pickle"
            reference result file
        checkpoint_1: str = "fine_filer_child_jobs.pickle"
            check point file that contain all children job objects
        work_dir: str = "./fine_filter"
            the working dir
        child_result_fname = "result.pickle"
            result fname relative to the children dir (i.e., f"{work_dir}/group_{i}/")
        mutant_fname = "mutants.pickle"
            mutant fname relative to the children dir (i.e., f"{work_dir}/group_{i}/")
        child_main_fname = "child_main.py"
            main script fname relative to the work dir
        kwargs_fname = "child_main_kwargs.pickle"
            fname of the kwargs for the main script relative to the work dir

    Return:
        the mutant space after screening."""
    if ncaa_param_lib_path is None:
        ncaa_param_lib_path = f"{work_dir}/ncaa_lib"
    num_mut_each_grp = mh.calc_average_task_num(len(mutants), shrapnel_groups)
    child_main_path = f"{work_dir}/{child_main_fname}"
    kwargs_file = f"{work_dir}/{kwargs_fname}"

    if not Path(checkpoint_1).exists():
        fs.safe_mkdir(work_dir)
        # 1. make main script
        # make kwargs
        child_main_kwargs = {
            "wt_stru" : stru,
            "atom_1" : stru.get(atom_1),
            "atom_2" : stru.get(atom_2),
            "ligand_chrg_spin_mapper" : ligand_chrg_spin_mapper,
            "ncaa_param_lib_path" : os.path.abspath(ncaa_param_lib_path),
            "cpu_job_config" : shrapnel_cpujob_config,
            "gpu_job_config" : shrapnel_gpujob_config,
            "md_constraints" : md_constraints,
            "md_length" : md_length,
            "md_parallel_runs" : md_parallel_runs,
            "ef_region_pattern" : ef_region_pattern,
        }
        save_obj(child_main_kwargs, kwargs_file)
        # make script
        save_func_to_main(_fine_filter_child_main, os.path.abspath(kwargs_file), child_main_path)
        # 2. seperate mutants in groups & make jobs
        child_jobs = []
        assigned = 0
        for grp_id, num_mut in enumerate(num_mut_each_grp):
            # make dir
            grp_path = f"{work_dir}/group_{grp_id}"
            fs.safe_mkdir(grp_path)
            # dist mutants
            grp_mut = mutants[assigned:assigned+num_mut]
            assigned += num_mut
            grp_mut_pickle = f"{grp_path}/{mutant_fname}"
            save_obj(grp_mut, grp_mut_pickle)
            # make job
            gpu_partition = _find_gpu_partition(grp_id, shrapnel_gpu_partition_mapper)
            sub_script_path = os.path.abspath(f"{grp_path}/submit_main.sh")
            child_jobs.append(
                _make_child_job(
                    grp_id = grp_id,
                    cluster_job_config = shrapnel_child_job_config,
                    child_main_path = os.path.abspath(child_main_path),
                    mutant_fname = mutant_fname,
                    child_result_fname = child_result_fname,
                    gpu_partition = gpu_partition,
                    sub_dir = grp_path,
                    sub_script_path = sub_script_path,
                )
            )

        # 3. submit jobs
        save_obj(child_jobs, checkpoint_1)
    else:
        child_jobs: List[ClusterJob] = load_obj(checkpoint_1)
        new_child_jobs = []
        # 1. analyze remaining child_jobs (recycle old one) may benefit from having a mimo
        for job in child_jobs:
            if not job.is_complete():
                new_child_jobs.append(job)
        child_jobs = new_child_jobs # this will include failing, pending, running jobs
        # 2. update the checkpoint
        save_obj(child_jobs, checkpoint_1)

    jobs_remain = ClusterJob.wait_to_array_end(child_jobs, shrapnel_check_period, shrapnel_child_array_size)

    if jobs_remain:
        _LOGGER.error("some children jobs didn't finish normally. they are:")
        _LOGGER.error("\n".join([j.sub_dir for j in jobs_remain]))
        raise ShrapnelChildError

    # 4. summarize result
    fine_metrics = {}
    for grp_id, num_mut in enumerate(num_mut_each_grp):
        grp_path = f"{work_dir}/group_{grp_id}"
        result_path = f"{grp_path}/{child_result_fname}"
        child_result = load_obj(result_path)
        fine_metrics.update(child_result)
    save_obj(fine_metrics, out_ref_path_1)

    # 5. screen based on fine metrics
    result_mutant_space = _screen_w_fine_metrics(mutants, fine_metrics, out_ref_path_2)

    return result_mutant_space

def _fine_filter_child_main(
        sys_argvs,
        wt_stru, # we dont do type hinting here as they will not be imported at the time of func def in the saved main script
        atom_1, atom_2,
        ligand_chrg_spin_mapper,
        ncaa_param_lib_path: str,
        cpu_job_config,
        gpu_job_config,
        md_constraints,
        md_length: float,
        md_parallel_runs: int,
        ef_region_pattern: str,
    ):
    """the child main script of the shrapnel treatment
    in fine_filter
    This will be the content of the main function with
    `-m mutant_file -p gpu_partition -o result` from cmdline."""
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
        if arg == "-m":
            mutants: List[List[Mutation]] = load_obj(sys_argvs[i+1])
        if arg == "-p":
            gpu_partition: str = sys_argvs[i+1]
        if arg == "-o":
            result_path: str = sys_argvs[i+1]
    # type hinting
    wt_stru: Structure
    ligand_chrg_spin_mapper: Dict
    atom_1: Atom
    atom_2: Atom
    cpu_job_config: Dict
    gpu_job_config: Dict
    md_constraints: List[Callable[[Structure], StructureConstraint]]
    # re-run
    result_dict = defaultdict(dict)
    if os.path.exists(result_path):
        result_dict: Dict = load_obj(result_path)

    # 1. prepare
    remove_hydrogens(wt_stru, polypeptide_only=True)
    protonate_stru(wt_stru, protonate_ligand=False)

    # 2. get ddg fold
    todo_mut = []
    for mut in mutants:
        mut_result_data: Dict = result_dict.get(tuple(mut), dict())
        if "ddg_fold" not in mut_result_data: # skip in re-run
            todo_mut.append(mut)
    if todo_mut:
        ddg_results = ddg_fold_of_mutants(
            wt_stru,
            todo_mut,
            num_iter = 10,
            cluster_job_config = cpu_job_config,
            relax_cluster_job_config = cpu_job_config,
        )
        for k, v in ddg_results.items():
            result_dict[k]["ddg_fold"] = v
        # save
        save_obj(result_dict, result_path)

    for i, mut in enumerate(mutants):
        tuple_mut = tuple(mut)
        mut_result_data: Dict = result_dict.get(tuple_mut, dict())
        # 3. mutate
        mutant_stru = mut_result_data.get("mutant_stru", None)
        if mutant_stru is None:
            mutant_dir = f"mutant_{i}" # as this will be executed under the grp_dir
            mutant_stru = mutate_stru(wt_stru, mut, engine="pymol")
            mutant_stru.assign_ncaa_chargespin(ligand_chrg_spin_mapper)
            remove_hydrogens(mutant_stru, polypeptide_only=True)
            protonate_stru(mutant_stru, protonate_ligand=False)
            result_dict[tuple_mut]["mutant_stru"] = mutant_stru
            # save
            save_obj(result_dict, result_path)
        
        # 4. MD
        trajs: List[StructureEnsemble] = mut_result_data.get("trajs", list())
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
            mut_constraints = []
            for cons in md_constraints:
                mut_constraints.append(cons(mutant_stru))
            new_trajs = equi_md_sampling(
                stru = mutant_stru,
                param_method = param_method,
                cluster_job_config = gpu_job_config,
                prod_constrain=mut_constraints,
                prod_time=md_length,
                record_period=md_length*0.01,
                work_dir=f"{mutant_dir}/MD/",
                parallel_runs=runs_left,
            )
            trajs.extend(new_trajs)
            # save
            result_dict[tuple_mut]["trajs"] = trajs
            save_obj(result_dict, result_path)
        
        mut_data = {
            "ef" : [],
            "spi" : [],
            "mmpbsa" : [],
            "mmgbsa" : [],
        }
        for replica_esm in trajs:
            replica_ef = []
            replica_spi = []
            for traj_stru in replica_esm.structures(remove_solvent=True):
                # 5. get dEF
                field_strength = ele_field_strength_at_along(
                    traj_stru, atom_1, atom_2, region_pattern=ef_region_pattern)
                replica_ef.append(field_strength)
                # 6. get SPI
                # spi = xxx()
                # replica_spi.append(spi)
            # 7. get MMPBSA
            # replica_mmpbsa, replica_mmgbsa = xxx()

            mut_data["ef"].append(replica_ef)
            # mut_data["spi"].append(replica_spi)
            # mut_data["mmpbsa"].append(replica_mmpbsa)
            # mut_data["mmgbsa"].append(replica_mmgbsa)

        mut_result_data.update(mut_data)
        result_dict[tuple_mut] = mut_result_data

    for mut in mutants:
        assert tuple(mut) in result_dict
        assert "ddg_fold" in result_dict[tuple(mut)]
        assert "ef" in result_dict[tuple(mut)]
        assert result_dict[tuple(mut)]["ef"]
        # assert "spi" in result_dict[tuple(mut)]
        # assert result_dict[tuple(mut)]["spi"]
        # assert "mmpbsa" in result_dict[tuple(mut)]
        # assert result_dict[tuple(mut)]["mmpbsa"]
        # assert "mmgbsa" in result_dict[tuple(mut)]
        # assert result_dict[tuple(mut)]["mmgbsa"]

    save_obj(result_dict, result_path)

def _make_child_job(
        grp_id: int,
        cluster_job_config: Dict,
        child_main_path: str,
        mutant_fname: str,
        child_result_fname: str,
        gpu_partition: str,
        sub_dir: str,
        sub_script_path: str,
        ) -> ClusterJob:
    """make the ClusterJob for children runs
    under the shrapnel framework"""
    cluster = cluster_job_config["cluster"]
    res_keywords = cluster_job_config["res_keywords"]

    cmd = f"python -u {child_main_path} -m {mutant_fname} -p {gpu_partition} -o {child_result_fname} > {child_main_path}.out 2>&1"
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

def _screen_w_fine_metrics(
        mutants: List[List[Mutation]],
        fine_metrics: Dict[Tuple[Mutation], Dict],
        out_ref_path_2: str) -> List[List[Mutation]]:
    """screen the mutant space with the fine metrics.
    and save a ref pickle file for all resulting mutants with
    fine metricss"""

def workflow_puo(single_mut_ddg_file_path: str = None):
    """the PuO workflow"""
    # target system
    atom_1 = "C.901.C64"
    atom_2 = "D.902.N7"
    chain_sync_list = [("A", "B")]
    chain_index_mapper = {
                "A" : 0,
                "B" : 450,
            }
    ncaa_chrgspin = {
        "FAD" : (0,1),
        "ACP" : (0,1),
    }
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "account" : "yang_lab_csb",
            "partition" : "production",
            'walltime' : '10:00:00',
        }
    }
    # file paths
    checkpoint_1 = "mutant_space_1.pickle"
    checkpoint_2 = "mutant_space_2.pickle"
    checkpoint_3 = "mutant_space_3.pickle"

    # 0. create Structure()
    target_pdb = glob.glob("*pdb")[0]
    stru = sp.get_structure(target_pdb)
    stru.assign_ncaa_chargespin(ncaa_chrgspin)

    # 1. optimizer
    if not Path(checkpoint_1).exists():
        mutant_space_1 = ef_optimzer(
            stru, atom_1, atom_2, 10.0,
            chain_sync_list, chain_index_mapper,
            limit_size=False)
        save_obj(mutant_space_1, checkpoint_1)
    else:
        mutant_space_1 = load_obj(checkpoint_1)

    # 2. coarse filter (mutant_space_2 are in groups)
    if not Path(checkpoint_2).exists():
        mutant_space_2 = coarse_filter(
            stru, mutant_space_1,
            atom_1, atom_2,
            cluster_job_config = cluster_job_config,
            ddg_fold_file=single_mut_ddg_file_path,)
        save_obj(mutant_space_2, checkpoint_2)
    else:
        mutant_space_2 = load_obj(checkpoint_2)

    # 3. fine filter
    if not Path(checkpoint_3).exists():
        shrapnel_child_job_config = {
            "cluster" : cluster_job_config["cluster"],
            "res_keywords" : cluster_job_config["res_keywords"] | {"walltime" : "10-00:00:00"}
        }
        shrapnel_cpujob_config = {
            "cluster" : cluster_job_config["cluster"],
            "res_keywords" : cluster_job_config["res_keywords"] | {"walltime" : "2-00:00:00"},
        }
        shrapnel_gpujob_config = {
            "cluster" : cluster_job_config["cluster"],
            "res_keywords" : cluster_job_config["res_keywords"] | {
                "account" : "csb_gpu_acc",
                "partition" : "pascal",
                "walltime" : "5-00:00:00",
            }
        }
        md_constraint = [
            partial(stru_cons.create_distance_constraint,
                "C.901.C64", "D.902.N7", 2.4),]
        mutant_space_2 = list(itertools.chain.from_iterable(mutant_space_2))
        mutant_space_3 = fine_filter(
            stru = stru, 
            mutants= mutant_space_2,
            atom_1 = atom_1,
            atom_2 = atom_2,
            ligand_chrg_spin_mapper = {"ACP" : (0,1), "FAD" : (0,1)},
            md_constraints=md_constraint,
            md_length=1.0, # TODO change this back after the test
            md_parallel_runs=1,
            ef_region_pattern="chain A+B+C+D+E+F and (not resi 901+902)",
            shrapnel_child_job_config = shrapnel_child_job_config,
            shrapnel_cpujob_config = shrapnel_cpujob_config,
            shrapnel_gpujob_config = shrapnel_gpujob_config,
            shrapnel_child_array_size = 2, # TODO change this after the test
            shrapnel_groups = 100,
            shrapnel_gpu_partition_mapper = {
                (0, 40) : "pascal",
                (41, 80) : "turing",
                (81, 99) : "a6000x4",
            },
        )
        save_obj(mutant_space_3, checkpoint_3)
    else:
        mutant_space_3 = load_obj(checkpoint_3)

    return mutant_space_3

def main():
    workflow_puo("single_mutation_ddg.pickle")

if __name__ == "__main__":
    main()
