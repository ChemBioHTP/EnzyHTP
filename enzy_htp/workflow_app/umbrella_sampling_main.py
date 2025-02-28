from functools import partial
import glob
import pickle
import sys

from enzy_htp.preparation import protonate_stru, remove_hydrogens
from enzy_htp.mutation import assign_mutant, mutate_stru
from enzy_htp.geometry import md_simulation
from enzy_htp.quantum import single_point
from enzy_htp.analysis import bond_dipole, ele_field_strength_at_along, ele_stab_energy_of_bond
from enzy_htp import interface
from enzy_htp._interface.amber_interface import AmberParameter
from enzy_htp._interface.handle_types import MolDynResult
import enzy_htp.structure.structure_constraint as stru_cons
from enzy_htp.structure import PDBParser
from enzy_htp.chemical.level_of_theory import QMLevelOfTheory
from enzy_htp.core.clusters.accre import Accre

# workflow config
# I/O path
wt_pdb_path = glob.glob("*amber.pdb")[0]
wt_inpcrd_path = glob.glob("*.inpcrd")[0]
wt_prmtop_path = glob.glob("*.prmtop")[0]
ligand_chrg_spin_mapper = {"LIGAND" : (0,1)}

force_constant = float(sys.argv[1])
target_sugar = int(sys.argv[2])
print(f"force_constant: {force_constant} target_sugar: {target_sugar}")
stru = PDBParser().get_structure(wt_pdb_path)
stru.assign_ncaa_chargespin(ligand_chrg_spin_mapper)
cons = [
    partial(
        stru_cons.create_group_distance_constraint,
        "resi 55 & n. C+CA+N", f"resi {target_sugar} & not elem H", topology=stru,
        params={
            "amber" : {
                "ialtd" : 0,
                "r1" : "x-1",
                "r2" : "x",
                "r3" : "x",
                "r4" : "x+1",
                "rk2": force_constant, "rk3": force_constant,
                "rs_filepath" : "{mdstep_dir}/0.rs",}}
    ),
]
windows = [
    [14.6],
    [15.6],
    [16.6],
    [17.6],
    [18.6],
    [19.6],
    [20.6],
    [21.6],
    [22.6],
    [23.6],
    [24.6],
]
prod_time = 20.0
prod_temperature = 300.0
record_period = 0.2

# HPC job resources
md_hpc_job_config = {
    "cluster" : Accre(),
    "res_keywords" : {
        "account" : "csb_gpu_acc",
        "partition" : "a6000x4",
        'walltime' : '3-00:00:00',
    }
}

# iterate over windows 
first = True
for i, window in enumerate(windows):
    window_constraints = [con(target_value=value) for con, value in zip(cons, window)]
    cluster_job_config = md_hpc_job_config
    if first:
        params = AmberParameter(wt_inpcrd_path, wt_prmtop_path, stru.ncaa_chrgspin_mapper) 
        first = False
    else:
        params = AmberParameter(last_frame_path, wt_prmtop_path, stru.ncaa_chrgspin_mapper) # use last rst as input

    min_step  = interface.amber.build_md_step(
        name="min_micro",
        minimize=True,
        length=20000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=window_constraints)

    # heat_step = interface.amber.build_md_step(
    #     name="heat_nvt",
    #     length=0.05, # ns
    #     cluster_job_config=cluster_job_config,
    #     core_type="gpu",
    #     temperature=[(0, 0), (0.05*0.9, prod_temperature), (-1, prod_temperature)],
    #     constrain=window_constraints)

    equi_step = interface.amber.build_md_step(
        name="equi_npt",
        length=prod_time * 0.01,
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=prod_temperature,
        constrain=window_constraints)

    prod_step = interface.amber.build_md_step(
        name="prod_npt",
        length=prod_time,
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        restart=True,
        if_report=True,
        temperature=prod_temperature,
        record_period=record_period,
        constrain=window_constraints)

    _, md_result = md_simulation(
        stru = stru,
        param_method = None,
        steps=[min_step,equi_step,prod_step],
        params_in = params,
        job_check_period=10,
        work_dir=f"./MD_w{i}/"
    )

    last_prod_result: MolDynResult = md_result[0][-1]
    last_frame_path = last_prod_result.last_frame_file

