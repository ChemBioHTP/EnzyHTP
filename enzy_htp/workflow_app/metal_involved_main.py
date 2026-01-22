from functools import partial
import pickle
from enzy_htp.preparation import protonate_stru, remove_hydrogens
from enzy_htp.mutation import assign_mutant, mutate_stru
from enzy_htp.geometry import equi_md_sampling
from enzy_htp import interface
import enzy_htp.structure.structure_constraint as stru_cons
from enzy_htp.structure import PDBParser
from enzy_htp.core.clusters.accre_r9 import AccreR9

# workflow config
# I/O path
wt_pdb_path = "3pzw_sub.pdb"
result_path = "result.pickle"

# define the charge spin for ligands and modified AAs
ligand_chrg_spin_mapper = {"FE11" : (3,6), "O11" : (-1,1), "LIG" :(-1,1), "HD1" : (0,1), "HD2" : (0,1), "HD3" : (0,1), "IE1" : (-1,1), "AN1" : (0,1)}
mod_aa_mainchain_atom_mapper = {"AN1": ["N", "CA", "C"],
                                "HD1": ["N", "CA", "C"],
                                "HD2": ["N", "CA", "C"],
                                "HD3": ["N", "CA", "C"],
                                "IE1": ["N", "CA", "C"]
                                } 

constraints=[]

# set up your ACCRE info
accre_res_account = "csb_gpu_acc"
accre_gpu_queue = "batch_gpu"

# ============================
# ==== main script lines =====
md_hpc_job_config = {
    "cluster" : AccreR9(),
    "res_keywords" : {
        "account" : accre_res_account, 
        "partition" : accre_gpu_queue,
        "node_cores" : "nvidia_rtx_a6000:1",
    }
}
result_dict = {}

# 1. create Structure()
wt_stru = PDBParser().get_structure(wt_pdb_path)

# 2. prepare
# remove_hydrogens(wt_stru, polypeptide_only=True)
# protonate_stru(wt_stru, protonate_ligand=True, )

# 3. create mutant library
mutant_pattern =  "WT" #, {L546A}, {L754A}, {I552A}, "{I553G}, {L546A,L754A}"
mutants = assign_mutant(wt_stru, mutant_pattern)

for i, mut in enumerate(mutants):
    mutant_result = []
    mutant_dir = f"mutant_ff14sb_{i}"

# 4. mutate Structure()
    mutant_stru = mutate_stru(wt_stru, mut, engine="pymol")
    mutant_stru.assign_ncaa_chargespin(ligand_chrg_spin_mapper)
    mutant_stru.assign_mod_aa_mainchain(mod_aa_mainchain_atom_mapper)

    # remove_hydrogens(mutant_stru, polypeptide_only=True)
    # protonate_stru(mutant_stru, protonate_ligand=True)
    # PDBParser().save_structure("protonated.pdb", mutant_stru)

# 5. sampling with MD
    param_method = interface.amber.build_md_parameterizer(
        ncaa_param_lib_path=f"ncaa_lib_14",
        force_fields=[
            # "leaprc.protein.ff19SB", # ff19SB is different from ff14SB by the CX XC type difference for CA
            "leaprc.protein.ff14SB",
            "leaprc.gaff",
            # "leaprc.water.opc",
            "leaprc.water.tip3p",
        ], 
        charge_method="AM1BCC",
        additional_tleap_lines=[
            'addAtomTypes {{ "M1"  "F" "sp3" }{ "Y1"  "N" "sp3" }{ "Y2"  "N" "sp3" }{ "Y3"  "N" "sp3" }{ "Y4"  "O" "sp3" }{ "Y5"  "O" "sp3" }{ "Y6"  "O" "sp3" }}',
            "loadamberparams frcmod.ionslm_126_opc",
            "bond a.499.NE2 a.840.FE1",
            "bond a.504.NE2 a.840.FE1",
            "bond a.690.NE2 a.840.FE1",
            "bond a.694.OD1 a.840.FE1",
            "bond a.839.OXT a.840.FE1",
            "bond a.840.FE1 a.841.O",
            "bond a.498.C a.499.N",
            "bond a.499.C a.500.N",
            "bond a.503.C a.504.N",
            "bond a.504.C a.505.N",
            "bond a.689.C a.690.N",
            "bond a.690.C a.691.N",
            "bond a.693.C a.694.N",
            "bond a.694.C a.695.N",
            "bond a.838.C a.839.N",
        ],
    )
    mut_constraints = [
        cons(topology=mutant_stru) for cons in constraints
    ]

    md_result = equi_md_sampling(
        stru = mutant_stru,
        param_method = param_method,
        cluster_job_config = md_hpc_job_config,
        job_check_period=10,
        prod_constrain=mut_constraints,
        prod_time= 100.0, #ns
        record_period= 0.1,
        work_dir=f"{mutant_dir}/MD/",
        #parallel_runs=num_rep,
        #prod_temperature=prod_temperature,
    )

    for replica_esm in md_result:
        replica_result = []

# save the result
with open(result_path, "wb") as of:
    pickle.dump(result_dict, of)
