"""Submodule contains code for charactization of thermostability of the enzyme.
+ ddg_fold_of_mutant()
    calculate the change of dg_fold upon mutations in a mutant

Author: QZ Shao <shaoqz@icloud.com>

Date: 2024-02-12
"""
from typing import List, Dict, Tuple

from enzy_htp.mutation_class import Mutation, is_mutant_wt
from enzy_htp import interface
from enzy_htp._interface.handle_types import ddGFoldEngine
from enzy_htp.structure import Structure
from enzy_htp.core import _LOGGER
from enzy_htp.core.job_manager import ClusterJob

def ddg_fold_of_mutants(
        stru: Structure,
        mutant_space: List[List[Mutation]],
        method: str = "rosetta_cartesian_ddg",
        parallel_method: str="cluster_job",
        cluster_job_config: Dict= None,
        job_check_period: int= 180, # s
        job_array_size: int= 100,
        work_dir: str="./ddG_fold",
        keep_in_file: bool=False,
        action_on_wt_at_start: bool = True,
        **kwargs,
        ) -> Dict[Tuple[Mutation], float]:
    """Calculate the change of dG_fold of the protein ({stru}) mutants in a mutant 
    space.

    Args:
        stru:
            the target molecule of the calculation represented as Structure()
        mutant_space:
            the target list of mutants to assess.
        method:
            the algorithm of the ddG_fold calculation. (see Details)
        parallel_method:
            the method to parallelize the multiple runs as well as executed
            a single run.
        cluster_job_config: 
            the config for cluster_job if it is used as the parallel method.
        job_check_period:
            the time cycle for update job state change (Unit: s)
        job_array_size:
            how many jobs are allowed to submit simultaneously. (0 means all -> len(inp))
            (e.g. 5 for 100 jobs means run 20 groups. All groups will be submitted and
            in each group, submit the next job only after the previous one finishes.)
        work_dir:
            the working dir that contains all the files in the process
        keep_in_file:
            whether keep the input file of the calculation
    Returns:
        A dictionary that record the ddG_fold of each mutant.
        {mutant : ddG_fold, ...}
    Details:
        The ddG_fold calculation is a popular way to assess the thermostability of a 
    mutant relative the wild-type, which has to be stable in the working temperature
    in most cases. It assesses the change of the free energy difference between the 
    folded state and the unfolded state of the protein.
    Available implemtations in this field:
        (supported as a "method" keyword)
        - Rosetta Cartesian ddG (https://www.rosettacommons.org/docs/latest/cartesian-ddG, https://pubs.acs.org/doi/10.1021/acs.jctc.6b00819)
            keyword: "rosetta_cartesian_ddg"
            0. relax WT strutcure
            1. repack (MC-based sampling) all mutating residues simultaneously.
            2. repack all residues that is directly in contact with the residues
               simultaneously.
            3. calculate Rosetta Energy difference between WT and these repacked mutants.
            4. repeat 1-3 N times and average the result
            * Backbone conformation change is not sampled
            * Rosetta Energy function is a semi-emperical energy function.
            * Fast and have acceptable accuracy.
        (not supported yet)
        - Free Energy Perturbation
            + slow
            + accurate
        - AI based method
            (under development)"""

    # get engine
    ddg_fold_engine = DDG_FOLD_ENGINE[method](
        cluster_job_config = cluster_job_config,
        work_dir = work_dir,
        keep_in_file = keep_in_file,
        **kwargs
    )

    # parallel methods
    if parallel_method == "cluster_job":
        if not cluster_job_config:
            _LOGGER.error("cluster_job is used but cluster_job_config is not given! "
                        "You need to at least specify the account and partition. "
                        "See test/geometry/test_sampling.py::test_equi_md_sampling_lv1() for an example.")
            raise ValueError
        result = _parallelize_ddg_fold_with_cluster_job(
            stru, mutant_space, ddg_fold_engine,
            job_check_period,
            job_array_size,
            action_on_wt_at_start = action_on_wt_at_start,
        )
    elif parallel_method is None:
        result = _serial_ddg_fold(
            stru, mutant_space, ddg_fold_engine,
            action_on_wt_at_start = action_on_wt_at_start,
        )
    else:
        _LOGGER.error(f"{parallel_method} is not a supported parallel_method")
        raise ValueError        

    return result

def _parallelize_ddg_fold_with_cluster_job(
        stru: Structure,
        mutant_space: List[List[Mutation]],
        engine: ddGFoldEngine,
        job_check_period: int,
        array_size: int,
        action_on_wt_at_start: bool,
        ) -> Dict[Tuple[Mutation], float]:
    """The parallelization method: cluster_job.
    This method will utilize ARMer@EnzyHTP and make each calculation a ClusterJob and
    parallalize them in a job array"""
    result = {}
    job_list = []
    result_eggs = []
    # 0. action on WT
    if action_on_wt_at_start:
        stru = engine.action_on_wt(stru)

    # 1. prep jobs
    for mutant in mutant_space:
        if is_mutant_wt(mutant): # treatment on WT
            result[tuple(mutant)] = 0.0
            continue
        job, egg = engine.make_job(stru, mutant)
        job_list.append(job)
        result_eggs.append((mutant, egg))
    
    # 2. wait to end
    ClusterJob.wait_to_array_end(
        job_list,
        period=job_check_period,
        array_size=array_size,
    )
    
    # 3. translate eggs
    for mutant, egg in result_eggs:
        result[tuple(mutant)] = engine.translate(egg)

    return result


def _serial_ddg_fold(
        stru: Structure,
        mutant_space: List[List[Mutation]],
        engine: ddGFoldEngine,
        action_on_wt_at_start: bool,
        ) -> Dict[Tuple[Mutation], float]:
    """The serial running method
    This method runs calculations in a serial manner locally."""
    result = []
    # 0. action on WT
    if action_on_wt_at_start:
        stru = engine.action_on_wt(stru)

    # 1. run jobs
    for mutant in mutant_space:
        output = engine.run(stru, mutant)
        result[tuple(mutant)] = output
    
    return result


DDG_FOLD_ENGINE: Dict[str, ddGFoldEngine] = {
    "rosetta_cartesian_ddg" : interface.rosetta.build_cartesian_ddg_engine,
}
