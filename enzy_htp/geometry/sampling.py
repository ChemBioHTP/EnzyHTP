"""Define functions for gemoetry sampling of Structure(). These functions sample gemoetries of
the target start structure on the energy surface defined by a certain energy function.

Science API:
    + md_simulation
    + equi_sampling()?

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-7-30
"""
from typing import List

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import InconsistentMDEngine
from enzy_htp.core import job_manager
from enzy_htp.core.mol_dyn_result import MolDynResult
from enzy_htp._interface.handle_types import MolDynStep, MolDynParameterizer
from enzy_htp.structure import Structure

def equi_md_sampling():
    """This science API performs a production run of molecular dynamics simulation with the
    system equilibrated by several short md simulations."""
    pass

# == general building blocks ==
def md_simulation(stru: Structure,
                  param_method: MolDynParameterizer,
                  steps: List[MolDynStep],
                  parallel_runs: int=1,
                  parallel_method: str="cluster_job",
                  work_dir: str="./MD") -> List[List[MolDynResult]]:
    """This science API perform Molecular Dynamics simulation.

    The starting structure {stru} will be parameterized by the {param_method} and
    simulated by sequential {steps} of {parallel_runs} numbers of parallel runs. The
    output of each last(?) step of the {steps} of each parallel runs will be returned as
    a list.

    Args:
        stru: the starting structure
        param_method: the Parameterizer() used for parameterization. This is a
            special step that covert enzy_htp.Structure() to the input format
            MolDynStep takes. Normally it will be topology, initial coordinate,
            and MM parameters etc.
        steps: a list of steps each is a MolDynStep() that defines a molecular
            dynamics step.
        parallel_runs: the number of desired parallel runs of the steps.
        parallel_method: the method to parallelize the multiple runs
        work_dir: the directory that contains all the MD files input/intermediate/output

    Return:
        a list of results that each results from the last step (unless specificed in the any step) of a parallel run

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
        results = _parallelize_md_steps_with_cluster_job(parallel_runs, work_dir, steps, params)
    ## sequential
    if parallel_method is None:
        results = _serial_md_steps(parallel_runs, work_dir, steps, params)

    return results

def _parallelize_md_steps_with_cluster_job(
        parallel_runs: int,
        work_dir: str,
        steps: List[MolDynStep],
        params,
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
            # 2. make output
            if step.if_report:
                result_egg_ele.append((step, output))
        if not result_egg_ele:
            result_egg_ele = [(step, output)] # default add last step if non is specified
        job_list = type(step).try_merge_jobs(job_list) # TODO this is a classmethod that try merging jobs of this type

        job_array.append(job_list)
        result_eggs.append(result_egg_ele)  # eggs are filenames that can be translated to give birth actual data

    job_manager.wait_to_2d_array_end(job_array) # TODO make this function that deal with job array that contain lists of serial jobs

    for rep_md_result in result_eggs:
        rep_result_list = []
        for step, output in rep_md_result:
            rep_result_list.append(step.translate(output))
        results.append(rep_result_list)

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

