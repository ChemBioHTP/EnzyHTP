"""The template for EnzyHTP 2.0 main script. It contains
A workflow that calculates the following properties for mutants:
- bond dipole
- electronic field strength
- electrostatic stablization energy (dG_ele)

Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-03-15"""
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

# workflow config
# I/O path
wt_pdb_path = "KE_07_R7_2_S.pdb"
result_path = "ke_test_result.pickle"
ligand_chrg_spin_mapper = {"H5J" : (0,1)} # define the charge spin for ligands and modified AAs

# HPC job resources
md_hpc_job_config = {
    "cluster" : Accre(),
    "res_keywords" : {
        "account" : "csb_gpu_acc",
        "partition" : "pascal"
    }
}
qm_hpc_job_config = {
    "cluster" : Accre(),
    "res_keywords" : {
        "account" : "yang_lab_csb",
        "partition" : "production",
        'walltime' : '1-00:00:00',
    }
}
result_dict = {}

# 1. create Structure()
wt_stru = PDBParser().get_structure(wt_pdb_path)

# 2. prepare
remove_hydrogens(wt_stru, polypeptide_only=True)
protonate_stru(wt_stru, protonate_ligand=False)

# 3. create mutant library
mutant_pattern = "WT, r:1[resi 254 around 4 and not resi 101: all not self]*10"
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
        force_fields=[
            "leaprc.protein.ff14SB",
            "leaprc.gaff",
            "leaprc.water.tip3p",
        ],
    )
    mut_constraints = [
        stru_cons.create_distance_constraint("B.254.H2", "A.101.OE2", 2.4, mutant_stru),
        stru_cons.create_angle_constraint("B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, mutant_stru),
    ]

    md_result = equi_md_sampling(
        stru = mutant_stru,
        param_method = param_method,
        cluster_job_config = md_hpc_job_config,
        job_check_period=10,
        prod_constrain=mut_constraints,
        prod_time= 0.4, #ns
        record_period=0.1,
        work_dir=f"{mutant_dir}/MD/"
    )

    for replica_esm in md_result:
        replica_result = []

# 6. electronic structure
        qm_results = single_point(
            stru=replica_esm,
            engine="gaussian",
            method=QMLevelOfTheory( basis_set="3-21G", method="hf" ),
            regions=["resi 101+254"],
            cluster_job_config=qm_hpc_job_config,
            job_check_period=60,
            job_array_size=20,
            work_dir=f"{mutant_dir}/QM_SPE/",
        )

# 7. analysis
        for ele_stru in qm_results:
            this_frame_stru = ele_stru.geometry.topology
            atom_1 = this_frame_stru.get("B.254.CAE")
            atom_2 = this_frame_stru.get("B.254.H2")

        # bond dipole
            dipole = bond_dipole(
                ele_stru, atom_1, atom_2,
                work_dir=f"{mutant_dir}/bond_dipole/"
            )
        # EF
            field_strength = ele_field_strength_at_along(
                this_frame_stru, atom_1, atom_2, region_pattern="chain A and (not resi 101)"
            )
        # dGele
            dg_ele = ele_stab_energy_of_bond(dipole[0], field_strength)

            replica_result.append((dg_ele, dipole, field_strength))

        mutant_result.append(replica_result)

    result_dict[tuple(mut)] = mutant_result

# save the result
    with open(result_path, "wb") as of:
        pickle.dump(result_dict, of)

