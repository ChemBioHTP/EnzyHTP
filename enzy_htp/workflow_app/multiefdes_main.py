"""the main script for the MultiMutDes-EF workflow.
Specialized for the optimization of PuO.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-11
"""
import glob
from pathlib import Path
import pickle
from typing import Any, Dict, List, Tuple

from enzy_htp.analysis import d_ele_field_upon_mutation_coarse, ddg_fold_of_mutants
from enzy_htp.core.clusters.accre import Accre
from enzy_htp.mutation import assign_mutant
from enzy_htp.mutation_class import Mutation, get_involved_mutation, generate_from_mutation_flag
from enzy_htp.structure import Structure, PDBParser, Atom
from enzy_htp.core import _LOGGER
import numpy as np

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
        result[tuple(mutant)].append(mutant_ef_strength)

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

def fine_filter() -> List[List[Mutation]]:
    """the fine filter applied on the 2nd mutant space
    and assess dEF, ddG_fold, MMPBSA, and SPI. Gives the final mutant space."""
    # NOTE MMPBSA is important for some close range

def save_obj(obj: Any, out_path: str):
    """save the mutant space to the out_path as a checkpoint
    .pickle file"""
    with open(out_path, "wb") as of:
        pickle.dump(obj, of)

def load_obj(in_path: str):
    """load the mutant space to the in_path as a checkpoint
    .pickle file"""
    with open(in_path, "rb") as f:
        obj = pickle.load(f)
    return obj

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

    # 2. coarse filter
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
        mutant_space_3 = fine_filter(mutant_space_2)
        save_obj(mutant_space_3, checkpoint_3)
    else:
        mutant_space_3 = load_obj(checkpoint_3)

    return mutant_space_3

def main():
    workflow_puo("single_mutation_ddg.pickle")

if __name__ == "__main__":
    main()
