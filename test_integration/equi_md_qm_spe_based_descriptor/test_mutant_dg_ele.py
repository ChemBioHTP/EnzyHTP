"""This file contain integration tests that calculates the change of electrostatic
stablization energy upon mutants for different systems.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-25"""
import pytest
import os
import pickle
from typing import List, Tuple, Dict, Callable
from functools import partial

from enzy_htp.preparation import protonate_stru, remove_hydrogens
from enzy_htp.mutation import assign_mutant, mutate_stru
from enzy_htp.geometry import equi_md_sampling
from enzy_htp.quantum import single_point
from enzy_htp.analysis import (
    bond_dipole, 
    ele_field_strength_at_along, 
    ele_stab_energy_of_bond)
from enzy_htp import interface
import enzy_htp.structure.structure_constraint as stru_cons
from enzy_htp.structure import (
    Structure,
    StructureConstraint,
    PDBParser
)
from enzy_htp.chemical.level_of_theory import QMLevelOfTheory
from enzy_htp.core.clusters.accre import Accre

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../../test/test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()

def workflow(
        wt_stru: Structure,
        mutant_pattern: str,
        ligand_chrg_spin_mapper: Dict,
        md_constraints: List[Callable[[Structure], StructureConstraint]],
        md_length: float,
        qm_region_pattern: str,
        qm_level_of_theory: QMLevelOfTheory,
        target_bond: Tuple[str],
        ef_region_pattern: str,
        result_path: str,
        chain_sync_list: List = None,
        chain_index_mapper: Dict = None,
    ):
    """the target workflow for tests in this file.
    This workflow calculates the dGele of each mutant
    and save it into a pickle file.
    Args:
        wt_pdb:
            the wild-type PDB
        mutant_pattern:
            the pattern of target mutants
        md_constraint:
            the constraint applyed during the sampling
        md_length:
            the length of the MD
        qm_region_pattern:
            the pattern of the target QM region
        target_bond:
            the target bond in (atom_key, atom_key)
        result_path:
            the path of the result file"""
    result_dict = {}
    # handles re-run
    if os.path.exists(result_path):
        with open(result_path, "rb") as f:
            result_dict = pickle.load(f)
    # init
    bond_p1, bond_p2 = target_bond
    
    # prepare
    remove_hydrogens(wt_stru, polypeptide_only=True)
    protonate_stru(wt_stru, protonate_ligand=False)

    # mutation
    mutants = assign_mutant(wt_stru, mutant_pattern, chain_sync_list, chain_index_mapper)
    for i, mut in enumerate(mutants):
        if tuple(mut) in result_dict:
            continue
        mutant_dir = f"{WORK_DIR}mutant_{i}"
        mutant_result = []
        mutant_stru = mutate_stru(wt_stru, mut, engine="pymol")
        mutant_stru.assign_ncaa_chargespin(ligand_chrg_spin_mapper)
        remove_hydrogens(mutant_stru, polypeptide_only=True)
        protonate_stru(mutant_stru, protonate_ligand=False)

    # sampling
        param_method = interface.amber.build_md_parameterizer(
            ncaa_param_lib_path=f"{DATA_DIR}ncaa_lib",
            force_fields=[
                "leaprc.protein.ff14SB",
                "leaprc.gaff",
                "leaprc.water.tip3p",
            ],
        )
        md_cluster_job_config = {
            "cluster" : Accre(),
            "res_keywords" : {
                "account" : "csb_gpu_acc",
                "partition" : "pascal"
            }
        }
        mut_constraints = []
        for cons in md_constraints:
            mut_constraints.append(cons(mutant_stru))

        md_result = equi_md_sampling(
            stru = mutant_stru,
            param_method = param_method,
            cluster_job_config = md_cluster_job_config,
            job_check_period=10,
            prod_constrain=mut_constraints,
            prod_time=md_length,
            record_period=md_length*0.1,
            work_dir=f"{mutant_dir}/MD/"
        )        

    # electronic structure
        for replica_esm in md_result:
            replica_result = []
            qm_cluster_job_config = {
                "cluster" : Accre(),
                "res_keywords" : {
                    "account" : "yang_lab_csb",
                    "partition" : "production",
                    'walltime' : '1-00:00:00',
                }}
            qm_results = single_point(
                stru=replica_esm,
                engine="gaussian",
                method=qm_level_of_theory,
                regions=[qm_region_pattern],
                cluster_job_config=qm_cluster_job_config,
                job_check_period=60,
                job_array_size=20,
                work_dir=f"{mutant_dir}/QM_SPE/",
            )

    # analysis
    # dGele
            for ele_stru in qm_results: # TODO make each API also handle list of them
                dipole = bond_dipole( # take around 1.5 h for 100 frames.
                    ele_stru, bond_p1, bond_p2,
                    work_dir=f"{mutant_dir}/bond_dipole/")
                field_strength = ele_field_strength_at_along(
                    ele_stru.geometry.topology, bond_p1, bond_p2, region_pattern=ef_region_pattern)
                dg_ele = ele_stab_energy_of_bond(dipole[0], field_strength)

                replica_result.append((dg_ele, dipole, field_strength))

            mutant_result.append(replica_result)

        result_dict[tuple(mut)] = mutant_result

        # update the pickle
        with open(result_path, "wb") as of:
            pickle.dump(result_dict, of)

@pytest.mark.accre
@pytest.mark.temp
def test_kemp_elimiase():
    """test the workflow on a kemp elimiase"""
    wt_stru = sp.get_structure(f"{STRU_DATA_DIR}/KE_07_R7_2_S.pdb")
    mutant_pattern = "WT, r:2[resi 254 around 4 and not resi 101: all not self]*2"
    md_constraint = [
        partial(stru_cons.create_distance_constraint,
            "B.254.H2", "A.101.OE2", 2.4),
        partial(stru_cons.create_angle_constraint,
            "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0),]
    qm_level_of_theory = QMLevelOfTheory(
        basis_set="3-21G",
        method="hf",        
    )
    target_bond = (
        wt_stru.ligands[0].find_atom_name("CAE"),
        wt_stru.ligands[0].find_atom_name("H2")
    )

    workflow(
        wt_stru = wt_stru,
        mutant_pattern = mutant_pattern,
        ligand_chrg_spin_mapper = {"H5J" : (0,1)},
        md_constraints = md_constraint,
        md_length = 0.1, #ns
        qm_region_pattern = "resi 101+254",
        qm_level_of_theory = qm_level_of_theory,
        target_bond = target_bond,
        ef_region_pattern = "chain A and (not resi 101)",
        result_path = f"{WORK_DIR}ke_test_result.pickle", 
    )
