"""Define functions for gemoetry sampling of Structure(). These functions sample gemoetries of
the target start structure on the energy surface defined by a certain energy function.

Science API:
    + md_simulation
    + equi_md_sampling()

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-7-30
"""
from typing import List, Dict, Tuple

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import InconsistentMDEngine
from enzy_htp.core import job_manager
from enzy_htp.structure import Structure, StructureEnsemble, structure_constraint
from enzy_htp._interface.handle_types import (
    MolDynStep,
    MolDynParameterizer,
    MolDynParameter,
    MolDynResult)

def equi_md_sampling(stru: Structure,
                     param_method: MolDynParameterizer,
                     parallel_runs: int= 3,
                     parallel_method: str= "cluster_job", # TODO prepare_only for just export files and cmd
                     work_dir: str="./MD",
                     # config for steps
                     prod_time: float= 50.0, # ns
                     prod_temperature: float = 300.0, #K
                     prod_constrain: structure_constraint.StructureConstraint= None,
                     record_period: float= 0.5, # ns
                     cluster_job_config: Dict= None,
                     cpu_equi_step: bool= False,
                     cpu_equi_job_config: Dict= None,
                     job_check_period: int=210, # s
                     ) -> List[StructureEnsemble]:
    """This science API performs a production run of molecular dynamics simulation with the
    system equilibrated by several short md simulations from the starting {stru}
    (Basically md_simulation() with preset steps)
    min (micro) -> heat (NVT) -> equi (NPT) -> prod (NPT)
    Args:
        stru: 
            the starting structure
        param_method: 
            the Parameterizer() used for parameterization. This determines the engine.
        parallel_runs: 
            the number of desired parallel runs of the steps.
        parallel_method: 
            the method to parallelize the multiple runs
        work_dir: 
            the directory that contains all the MD files input/intermediate/output
        prod_time: 
            the simulation time in production step (unit: ns)
        prod_temperature: 
            the production temperature
        prod_constrain: 
            the constrain applied in the production step
        record_period: 
            the simulation time period for recording the geom. (unit: ns)
        cluster_job_config: 
            the config for cluster_job if it is used as the parallel method.
        cpu_equi_step: 
            whether use cpu for equi step
        cpu_equi_job_config: 
            the job config for the cpu equi step if specified
        job_check_period:
            the check period for wait_to_2d_array_end. Used when parallel_method='cluster_job'.
            (Unit: s, default: 210s)
    Returns:
        a list trajectories for each replica in StructureEnsemble format."""
    result = []
    # san check
    if parallel_method == "cluster_job":
        if not cluster_job_config:
            _LOGGER.error("cluster_job is used but cluster_job_config is not given! "
                          "You need to at least specify the account and partition. "
                          "See test/geometry/test_sampling.py::test_equi_md_sampling_lv1() for an example.")
            raise ValueError

    # 1. build steps
    parent_interface = param_method.parent_interface

    # 1.1 equi core
    equi_core = "gpu"
    equi_job_config = cluster_job_config
    if cpu_equi_step:
        equi_core = "cpu"
        equi_job_config = cpu_equi_job_config
        if not cpu_equi_job_config:
            _LOGGER.error("cpu_equi_step is used but cpu_equi_job_config is not given! "
                          "You need to at least specify the account and partition. ")
            raise ValueError

    freeze_backbone = structure_constraint.create_backbone_freeze(stru)
    min_step  = parent_interface.build_md_step(
        name="min_micro",
        minimize=True,
        length=20000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=[freeze_backbone, prod_constrain])

    heat_step = parent_interface.build_md_step(
        name="heat_nvt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=[(0, 0), (0.05*0.9, prod_temperature), (-1, prod_temperature)],
        constrain=[freeze_backbone, prod_constrain])

    equi_step = parent_interface.build_md_step(
        name="equi_npt",
        length=prod_time * 0.01,
        cluster_job_config=equi_job_config,
        core_type=equi_core,
        temperature=prod_temperature,
        constrain=[freeze_backbone, prod_constrain])

    equi_step = parent_interface.build_md_step(
        name="equi_npt_free_bb",
        length=prod_time * 0.01,
        cluster_job_config=equi_job_config,
        core_type=equi_core,
        temperature=prod_temperature,
        constrain=[prod_constrain])

    prod_step = parent_interface.build_md_step(
        name="prod_npt",
        length=prod_time,
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        restart=True,
        if_report=True,
        temperature=prod_temperature,
        record_period=record_period,
        constrain=[prod_constrain])

    # 2. run simulation
    params, md_result = md_simulation(
        stru, param_method,
        steps=[min_step, heat_step, equi_step, prod_step],
        parallel_runs=parallel_runs,
        parallel_method=parallel_method,
        work_dir=work_dir,
        job_check_period=job_check_period,)

    # 3. format output
    for rep_result in md_result:
        prod_result = rep_result[-1]
        result.append(StructureEnsemble(
            topology=params.topology_file,
            top_parser=params.topology_parser,
            coordinate_list=prod_result.traj_file,
            coord_parser=prod_result.traj_parser,))

    return result

# == general building blocks ==
def md_simulation(stru: Structure,
                  param_method: MolDynParameterizer,
                  steps: List[MolDynStep],
                  parallel_runs: int=1,
                  parallel_method: str="cluster_job",
                  work_dir: str="./MD",
                  job_check_period: int=210,) -> Tuple[MolDynParameter, List[List[MolDynResult]]]:
    """This science API perform Molecular Dynamics simulation.

    The starting structure {stru} will be parameterized by the {param_method} and
    simulated by sequential {steps} of {parallel_runs} numbers of parallel runs. The
    output of each last(?) step of the {steps} of each parallel runs will be returned as
    a list.

    Args:
        stru:
            the starting structure
        param_method:
            the Parameterizer() used for parameterization. This is a
            special step that covert enzy_htp.Structure() to the input format
            MolDynStep takes. Normally it will be topology, initial coordinate,
            and MM parameters etc.
        steps:
            a list of steps each is a MolDynStep() that defines a molecular
            dynamics step.
        parallel_runs:
            the number of desired parallel runs of the steps.
        parallel_method:
            the method to parallelize the multiple runs
        work_dir:
            the directory that contains all the MD files input/intermediate/output
        job_check_period:
            the check period for wait_to_2d_array_end. Used when parallel_method='cluster_job'.
            (Unit: s, default: 210s)

    Return:
        Tuple[
            the parameter object,
            a list of results of each step from each parallel run
        ]
        example: (params, [[MolDynResult, MolDynResult, ...], ...])

    Details:
        Molecular Dynamics (MD) simulate the motion of enzymes based on Newton's equation of motion.
    The forces between every interacting atoms are determined by molecular mechanical force fields. (or
    QM in QM/MM MD) Upon assigning the initial velocity to the initial structure (atomic position), the
    equation derives the new velocity and new atomic position after each time step and yields a trajectory
    , that is, a collection of positions at each time point, after a certain length of the simulation.
    There are many detailed problems in MD that won't be discussed here. (revent reviews:
    https://www.tandfonline.com/doi/full/10.2147/AABC.S70333, https://www.nature.com/articles/nsb0902-646)

    Common implemtations of MD for enzyme in the field are:
    (TODO do more search on this, complete the note for each packages)
        - Amber (https://ambermd.org/, https://ambermd.org/2012_wires.pdf, https://github.com/Amber-MD)
            Large community, may be the best GPU speed for MD.
            Contain QUICK that support GPU-accelerated QM for QM/MM MD.
            Not Free.
        - OpenMM (https://github.com/openmm/openmm, https://github.com/molmod/openmm-tutorial-msbs/tree/main)
            A Python library for MD. Python friendly. Colabed with Tinker and have a AMOEBA plugin.
            Doesn't support FF19SB in normal way. Can only convert from prmtop generated using ff19SB.
            Free and open-source
        - Tinker (https://dasher.wustl.edu/tinker/)
            Support GPU-accelerate AMOEBA polarized FF.
            Free.
        - CHARMM (https://www.charmm.org/)
            Support multisite lambda dynamics (MSLD)
            Not Free.
        - NAMD (https://www.ks.uiuc.edu/Research/namd/)
            Support BFEE2 for automated free energy calculation.
            Scale to 100,000+ cores. (according to https://youtu.be/_TiQYNWJwYg?list=PLGL4XGw8noUwrh16gsC9H_D03fED3IcHo&t=2229)
            Free.
        - Gromacs (https://www.gromacs.org/)
            Free and open-source.
        - M-Chem (https://pubmed.ncbi.nlm.nih.gov/37470065/)
            Not Free.
            Support good general force field for ligands
            Support better implemtation of polarizable force field
    """
    supported_parallel_method = ["cluster_job"]
    # I. san check
    #   - MD engine consistency
    for i, step in enumerate(steps):
        if step.engine != param_method.engine:
            _LOGGER.error(
                f"The engine of step #{i} ({step.engine}) does not match the parameterizer ({param_method.engine})!")
            raise InconsistentMDEngine
    #   - Parallel method suppoort
    if parallel_method and (parallel_method not in supported_parallel_method):
        _LOGGER.error(
            f"Parallel method: {parallel_method} is not in the supported list: {supported_parallel_method}")
        raise ValueError

    # II. make work dir
    fs.safe_mkdir(work_dir)

    # III. parameterize
    params = param_method.run(stru)

    # IV. run MD steps
    ## parallelize
    if parallel_method == "cluster_job":
        results = _parallelize_md_steps_with_cluster_job(parallel_runs, work_dir, steps, params, job_check_period)
    ## sequential
    if parallel_method is None:
        results = _serial_md_steps(parallel_runs, work_dir, steps, params)

    return params, results

def _parallelize_md_steps_with_cluster_job(
        parallel_runs: int,
        work_dir: str,
        steps: List[MolDynStep],
        params: MolDynParameter,
        period: int,
        ):
    """The MD parallelization method: cluster_job. (only used in md_simulation())
    This method will utilize ARMer@EnzyHTP and make each MD steps a ClusterJob and
    parallalize them in a 2d job array (since there are some sequential dependency in MD steps)"""
    job_array = []
    result_eggs = []
    results = []
    for i in range(parallel_runs):
        job_list = []
        result_egg_ele = []
        # create job path
        sub_work_dir = f"{work_dir}/rep_{i}"
        fs.safe_mkdir(sub_work_dir)
        output = None  # the output place holder; the output between steps are very different for different packages so it will prob also becomes a class

        for step in steps:
            # 1. make job list
            step.work_dir = sub_work_dir
            if output: # steps after will use output from the previous one
                job, output = step.make_job(output)
            else: # the 1st step
                job, output = step.make_job(params)
            job_list.append(job)
            # 2. make output (we need to translate all since error checking and cleaning is needed)
            result_egg_ele.append((step, output))

        job_list = type(step).try_merge_jobs(job_list)

        job_array.append(job_list)
        result_eggs.append(result_egg_ele)  # eggs are filenames that can be translated to give birth actual data

    job_manager.ClusterJob.wait_to_2d_array_end(job_array, period=period)

    for rep_md_result in result_eggs:
        rep_result_list = []
        for step, output in rep_md_result:
            rep_result_list.append(step.translate(output))
        results.append(rep_result_list)

    # clean up
    job_temp_files = set()
    for job_list in job_array:
        for job in job_list:
            job: job_manager.ClusterJob
            job_temp_files.add(job.sub_script_path)
            # job_temp_files.add(job.job_cluster_log) probably not.
    fs.clean_temp_file_n_dir(list(job_temp_files))

    return results

def _serial_md_steps(
        parallel_runs: int,
        work_dir: str,
        steps: List[MolDynStep],
        params,
        ):
    """The MD serial running method (only used in md_simulation())
    This method runs MD steps in a serial manner locally."""
    results = []
    for i in range(parallel_runs):
        # create job path
        sub_work_dir = f"{work_dir}/rep_{i}"
        fs.safe_mkdir(sub_work_dir)
        output = None
        result_ele = []
        for step in steps:
            step.work_dir = sub_work_dir
            if not output:
                output = step.run(params)
            else:
                output = step.run(output)
            if step.if_report:
                result_ele.append((step, output))
        if not result_ele:
            result_ele = [(step, output)] # default add last step if non is specified

        results.append([step.translate(output) for step, output in result_ele])

    return results

