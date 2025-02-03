"""The template for EnzyHTP 2.0 main script. It run MD simulations with given setting.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2025-01-14"""
from functools import partial
import pickle
from enzy_htp.preparation import protonate_stru, remove_hydrogens
from enzy_htp.mutation import assign_mutant, mutate_stru
from enzy_htp.geometry import equi_md_sampling
from enzy_htp.quantum import single_point
from enzy_htp.analysis import bond_dipole, ele_field_strength_at_along, ele_stab_energy_of_bond
from enzy_htp import interface
import enzy_htp.structure.structure_constraint as stru_cons
from enzy_htp.structure import PDBParser
from enzy_htp.chemical.level_of_theory import QMLevelOfTheory
from enzy_htp.core.clusters.accre import Accre

# ==============
# == settings ==
# target structure (currently you need to dock it before hand)
wt_pdb_path = "KE_07_R7_2_S.pdb"
# define the charge spin for ligands and modified AAs
ligand_chrg_spin_mapper = {
    "H5J" : (0,1)
}
# protonation of ligand could be funky. Turn this off and do it manually when having problem.
protonate_ligand = True 
# mutation (learn syntax from https://enzyhtp-doc.readthedocs.io/en/latest/sci_api_tutorial/assign_mutant.html#mutant-pattern)
mutant_pattern = "WT, r:2[resi 254 around 4 and not resi 101: all not self]*2" # leave just "WT" if no mutations are wanted
# force field
force_fields = [
    "leaprc.protein.ff14SB",
    "leaprc.gaff",
    "leaprc.water.tip3p",
]
# constraints
constraints = [
    partial(stru_cons.create_distance_constraint, "B.254.H2", "A.101.OE2", 2.4),
    partial(stru_cons.create_angle_constraint, "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0),
]
# length on MD simulation
prod_time = 10.0 # unit: ns
# number of replica
num_rep = 3 
# temp of MD
prod_temperature = 300.0 #unit: K
# set up your ACCRE info
accre_res_account = "csb_gpu_acc"
accre_gpu_queue = "pascal"


# ============================
# ==== main script lines =====
md_hpc_job_config = {
    "cluster" : Accre(),
    "res_keywords" : {
        "account" : accre_res_account, 
        "partition" : accre_gpu_queue
    }
}

# 1. create Structure()
wt_stru = PDBParser().get_structure(wt_pdb_path)

# 2. prepare
remove_hydrogens(wt_stru, polypeptide_only=True)
protonate_stru(wt_stru, protonate_ligand=protonate_ligand)

# 3. create mutant library

mutants = assign_mutant(wt_stru, mutant_pattern)

for i, mut in enumerate(mutants):
    mutant_result = []
    mutant_dir = f"mutant_{i}"

# 4. mutate Structure()
    mutant_stru = mutate_stru(wt_stru, mut, engine="pymol")
    mutant_stru.assign_ncaa_chargespin(ligand_chrg_spin_mapper)

    remove_hydrogens(mutant_stru, polypeptide_only=True)
    protonate_stru(mutant_stru, protonate_ligand=False)

# 5. sampling with MD
    param_method = interface.amber.build_md_parameterizer(
        ncaa_param_lib_path=f"ncaa_lib",
        force_fields=force_fields,
    )
    mut_constraints = [
        cons(topology=mutant_stru) for cons in constraints
    ]

    md_result = equi_md_sampling(
        stru = mutant_stru,
        param_method = param_method,
        cluster_job_config = md_hpc_job_config,
        job_check_period=30,
        prod_constrain=mut_constraints,
        prod_time= prod_time, #ns
        work_dir=f"{mutant_dir}/MD/",
        parallel_runs=num_rep,
        prod_temperature=prod_temperature,
    )


