"""MVP script for doing umbrella sampling in parallel in EnzyHTP."""
from functools import partial
import glob
import pickle
import sys
import multiprocessing

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
from enzy_htp import config as eh_config
from enzy_htp.chemical.level_of_theory import QMLevelOfTheory
from enzy_htp.core.clusters.accre import Accre

def sequential_umbrella_sampling():
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
        [14.6], # 1st attempt
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
    first = True # 1st attempt

    for i, window in enumerate(windows):
        i = i # 1st attempt

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

def sequential_umbrella_sampling_w_start_point():
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
        [25.6],
        [26.6],
        [27.6],
        [28.6],
        [29.6],
        [30.6],
        [31.6],
        [32.6],
        [33.6],
        [34.6],
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
    first = False # 2nd attempt
    last_frame_path = "MD_w10/rep_0/prod_npt.rst" # 2nd attempt

    for i, window in enumerate(windows):
        i = i+11 # 2nd attempt

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

def parallel_umbrella_sampling_w_start_points():
    # workflow config
    # I/O path
    wt_pdb_path = glob.glob("*amber.pdb")[0]
    wt_prmtop_path = glob.glob("*.prmtop")[0]
    ligand_chrg_spin_mapper = {"LIGAND" : (0,1)}

    force_constant = float(sys.argv[1])
    target_sugar = int(sys.argv[2])
    partition = "a6000x4"
    if len(sys.argv) >= 4:
        partition = str(sys.argv[3])
        
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
    window_start_point = [
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
        [25.6],
        [26.6],
        [27.6],
        [28.6],
        [29.6],
        [30.6],
        [31.6],
        [32.6],
        [33.6],
        [34.6],
    ]
    prod_time = 20.0
    prod_temperature = 300.0
    record_period = 0.2

    # HPC job resources
    md_hpc_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "account" : "csb_gpu_acc",
            "partition" : partition,
            'walltime' : '3-00:00:00',
        }
    }
    task_args = [
        (i, start_point, cons, wt_prmtop_path, stru, md_hpc_job_config, prod_time, prod_temperature, record_period)
        for i, start_point in enumerate(window_start_point)
    ]
    with multiprocessing.Pool(processes=21) as pool:
        pool.map(grow_each_window, task_args)

def grow_each_window(args):
    i, start_point, cons, wt_prmtop_path, stru, md_hpc_job_config, prod_time, prod_temperature, record_period = args
    
    eh_config.system.SCRATCH_DIR = f"./scratch_{i}/"
    start_point_dir = f"MD_w{i}"
    last_frame_path = glob.glob(f"{start_point_dir}/rep_0*/prod_npt.rst")[0]
    print(f"starting from {last_frame_path}")
    filling_windows = [[start_point[0]+(j+1)*0.25] for j in range(3)]

    for k, window in enumerate(filling_windows):
        params = AmberParameter(last_frame_path, wt_prmtop_path, stru.ncaa_chrgspin_mapper) # use last rst as input
        window_constraints = [con(target_value=value) for con, value in zip(cons, window)]
        cluster_job_config = md_hpc_job_config

        min_step  = interface.amber.build_md_step(
            name="min_micro",
            minimize=True,
            length=20000, # cycle
            cluster_job_config=cluster_job_config,
            core_type="gpu",
            constrain=window_constraints)

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
            work_dir=f"./MD_w{i}_{k}/"
        )

        last_prod_result: MolDynResult = md_result[0][-1]
        last_frame_path = last_prod_result.last_frame_file

def a_pymol_bug_related_to_this_app():
    """see https://github.com/schrodinger/pymol-open-source/issues/436"""
    import multiprocessing
    import pymol2

    def grow_each_window(args):
        print("start")
        call_read_pdbstr()
        print("pass")

    def call_read_pdbstr():
        pms = pymol2.PyMOL()
        pms.start()
        pdb_str = "ATOM      1  N   ASP A   1      45.117  64.639  61.562  1.00  0.00"
        pms.cmd.read_pdbstr(pdb_str, "obj1")
        pms.stop()

    if __name__ == "__main__":
        # multiprocessing.set_start_method("spawn")
        call_read_pdbstr()
        with multiprocessing.Pool(processes=4) as pool:
            pool.map(grow_each_window, range(3))

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    parallel_umbrella_sampling_w_start_points()
