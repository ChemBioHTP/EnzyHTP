"""Submodule contains code for calculations the binding energy of a ligand in the enzyme.
The energy can be either enthalpy or Gibbs free energy depending on the algorithm used
+ binding_energy()
    calculate the binding energy of a ligand in a protein from a Structure or an ensemble
    of Structure()s

Author: QZ Shao <shaoqz@icloud.com>

Date: 2024-05-30
"""
from typing import List, Dict, Tuple, Union, Callable
from plum import dispatch

from enzy_htp import interface
from enzy_htp.structure.structure_selection import select_stru
from enzy_htp.structure import Structure, Ligand, StructureEnsemble
from enzy_htp.core import _LOGGER
from enzy_htp.core.job_manager import ClusterJob

def binding_energy(
        # Input Data
        stru: Union[Structure, StructureEnsemble],
        ligand: str,
        # Config
        method: str = "mmpbgbsa_amber",
        cluster_job_config: Dict= None,
        job_check_period: int= 30, # s
        work_dir: str="./binding",
        keep_in_file: bool=False,
        **kwargs,
        ) -> Union[float, List[float]]:
    """Calculate the binding energy of {ligand} in {stru}

    Args:
        stru:
            the target protein as the receptor of the calculation represented as Structure()
        ligand:
            the target ligand of the calculation represented as a selection pattern.
            Note that the ligand has to be part of Structure().
            Note that the ligand can be a small molecule or a protein   
        method:
            the algorithm of the binding energy calculation. (see Details)
        cluster_job_config: 
            the config for cluster_job. (also enables running via a ClusterJob/ARMer)
        job_check_period:
            the time cycle for update job state change (Unit: s)
        work_dir:
            the working dir that contains all the files in the process
        keep_in_file:
            whether keep the input file of the calculation
    (when method="mmpbsa")
        igb:
        TODO

    Returns:
        A float or a list of float (depending on whether an ensemble is given) for the value of the binding energy

    Details:
        The binding energy assesses how favorable would a ligand bind into a protein. It is the energy 
        difference between the bound state (protein-ligand) and the unbound state.

    Available implemtations in this field:
        (supported as a "method" keyword)
        - MMPB/GBSA-1A (https://www.tandfonline.com/doi/full/10.1517/17460441.2015.1032936, https://ambermd.org/tutorials/advanced/tutorial3/index.php, https://pubs.acs.org/doi/10.1021/ar000033j)
            keyword: "mmpbgbsa_amber"
            calculate the MMPB/GBSA binding energy from a single trajectory/StructureEnsemble of the complex using AmberTools
            For each frame from the trajectory, 
            MM:
                the MM energy of the dry complex, dry ligand, and dry protein is calculated 
                and the dry binding energy is calculated by E_dry = E_solvated_complex - E_dry_complex - E_dry_ligand
            PB/GB:
                implicit solvent model is used to calculated the solvation energy for each component
            SA:
                the surface area approximation is used for calculating the hydrophobic component of the solvation free energy
            The total binding energy in solvent is E_dry + E_solv_complex - E_solv_ligand - E_solv_protein
            * Most MMPB/GBSA calculation result approximates the binding enthalpy
            * MMPB/GBSA-1A is bad at calculeting protein-protein interaction (https://pubs.rsc.org/en/content/articlelanding/2016/cp/c6cp03670h)

        (not supported yet)
        - Rosetta Binding Energy
        - Free Energy Perturbation
            + slow
            + accurate
        - AI based method
            (under development)"""

    # # get engine
    # ddg_fold_engine = DDG_FOLD_ENGINE[method](
    #     cluster_job_config = cluster_job_config,
    #     work_dir = work_dir,
    #     keep_in_file = keep_in_file,
    #     **kwargs
    # )

    # # parallel methods
    # if parallel_method == "cluster_job":
    #     if not cluster_job_config:
    #         _LOGGER.error("cluster_job is used but cluster_job_config is not given! "
    #                     "You need to at least specify the account and partition. "
    #                     "See test/geometry/test_sampling.py::test_equi_md_sampling_lv1() for an example.")
    #         raise ValueError
    #     result = _parallelize_ddg_fold_with_cluster_job(
    #         stru, mutant_space, ddg_fold_engine,
    #         job_check_period,
    #         job_array_size,
    #         action_on_wt_at_start = action_on_wt_at_start,
    #     )
    # elif parallel_method is None:
    #     result = _serial_ddg_fold(
    #         stru, mutant_space, ddg_fold_engine,
    #         action_on_wt_at_start = action_on_wt_at_start,
    #     )
    # else:
    #     _LOGGER.error(f"{parallel_method} is not a supported parallel_method")
    #     raise ValueError        

    return result

# def _parallelize_ddg_fold_with_cluster_job(
#         stru: Structure,
#         mutant_space: List[List[Mutation]],
#         engine: ddGFoldEngine,
#         job_check_period: int,
#         array_size: int,
#         action_on_wt_at_start: bool,
#         ) -> Dict[Tuple[Mutation], float]:
#     """The parallelization method: cluster_job.
#     This method will utilize ARMer@EnzyHTP and make each calculation a ClusterJob and
#     parallalize them in a job array"""
#     result = {}
#     job_list = []
#     result_eggs = []
#     # 0. action on WT
#     if action_on_wt_at_start:
#         stru = engine.action_on_wt(stru)

#     # 1. prep jobs
#     for mutant in mutant_space:
#         if is_mutant_wt(mutant): # treatment on WT
#             result[tuple(mutant)] = 0.0
#             continue
#         job, egg = engine.make_job(stru, mutant)
#         job_list.append(job)
#         result_eggs.append((mutant, egg))
    
#     # 2. wait to end
#     ClusterJob.wait_to_array_end(
#         job_list,
#         period=job_check_period,
#         array_size=array_size,
#     )
    
#     # 3. translate eggs
#     for mutant, egg in result_eggs:
#         result[tuple(mutant)] = engine.translate(egg)

#     return result


# def _serial_ddg_fold(
#         stru: Structure,
#         mutant_space: List[List[Mutation]],
#         engine: ddGFoldEngine,
#         action_on_wt_at_start: bool,
#         ) -> Dict[Tuple[Mutation], float]:
#     """The serial running method
#     This method runs calculations in a serial manner locally."""
#     result = []
#     # 0. action on WT
#     if action_on_wt_at_start:
#         stru = engine.action_on_wt(stru)

#     # 1. run jobs
#     for mutant in mutant_space:
#         output = engine.run(stru, mutant)
#         result[tuple(mutant)] = output
    
#     return result


BINDING_ENERGY_ENGINE: Dict[str, Callable] = {
    "mmpbgbsa_amber" : interface.amber.get_mmpbgbsa_energy,
}
