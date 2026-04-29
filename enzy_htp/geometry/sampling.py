"""Define functions for gemoetry sampling of Structure(). These functions sample gemoetries of
the target start structure on the energy surface defined by a certain energy function.

Science API:
    + md_simulation()
    + equi_md_sampling()
    + deployable_equi_md_sampling()
    + umbrella_sampling()

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-7-30
"""
import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, List, Dict, Optional, Tuple, Union

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import InconsistentMDEngine
from enzy_htp.core import job_manager
from enzy_htp.structure import Structure, StructureEnsemble
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp._interface.handle_types import (
    MolDynStep,
    MolDynParameterizer,
    MolDynParameter,
    MolDynResult)


# ---------------------------------------------------------------------------
# Umbrella sampling result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class UmbrellaWindowResult:
    """Result produced by a single umbrella sampling window.

    Attributes:
        window_index:   Zero-based index of this window within the CV target list.
        window_target:  The RC target value for this window (in ``cv.unit``).
        md_results:     Raw MD output: ``List[List[MolDynResult]]`` where the
                        outer list spans parallel replicas and the inner list
                        spans MD steps.
        metadata:       Arbitrary key→value metadata; always includes:
                        ``cv`` (``cv.to_dict()``), ``engine_payload`` (the
                        dict returned by ``cv.serialize_for_engine``), and
                        ``work_dir``.
    """
    window_index: int
    window_target: float
    md_results: List[List[MolDynResult]]
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class UmbrellaSamplingResult:
    """Aggregated result from :func:`umbrella_sampling`.

    Attributes:
        windows:         List of per-window results in target order.
        cv_dict:         ``cv.to_dict()`` snapshot captured at call time.
        window_targets:  The ordered list of RC target values that were run.
        work_dir:        Root working directory used for the run.
    """
    windows: List[UmbrellaWindowResult]
    cv_dict: Dict[str, Any]
    window_targets: List[float]
    work_dir: str

def equi_md_sampling(stru: Structure,
                     param_method: MolDynParameterizer, # TODO support using engine + kwarg to specify
                     parallel_runs: int= 3,
                     parallel_method: str= "cluster_job",
                     work_dir: str="./MD",
                     # config for steps
                     prod_time: float= 50.0, # ns
                     prod_temperature: float = 300.0, #K
                     prod_constrain: List[stru_cons.StructureConstraint]= None,
                     record_period: float= 0.5, # ns
                     cluster_job_config: Dict= None,
                     cpu_equi_step: bool= False,
                     cpu_equi_job_config: Dict= None,
                     dont_freeze_bb_in_min: bool= False,
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
        dont_freeze_bb_in_min:
            the option allows one to remove backbone freeze during minimization
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

    param_method, (
        min_step, heat_step, equi_step_1, equi_step_2, prod_step
    ) = _process_equi_md_sampling_arguments(
        stru = stru,
        param_method = param_method,
        prod_time = prod_time,
        prod_temperature = prod_temperature,
        prod_constrain = prod_constrain,
        record_period = record_period,
        cluster_job_config = cluster_job_config,
        cpu_equi_step = cpu_equi_step,
        cpu_equi_job_config = cpu_equi_job_config,
        dont_freeze_bb_in_min = dont_freeze_bb_in_min,
    )

    # 2. run simulation
    params, md_result = md_simulation(
        stru, param_method,
        steps=[min_step, heat_step, equi_step_1, equi_step_2, prod_step],
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

def deployable_equi_md_sampling(
        stru: Structure,
        param_method: MolDynParameterizer,
        parallel_runs: int= 3,
        work_dir: str="./MD",
        # config for steps
        prod_time: float= 50.0, # ns
        prod_temperature: float = 300.0, #K
        prod_constrain: List[stru_cons.StructureConstraint]= None,
        record_period: float= 0.5, # ns
        cluster_job_config: Dict= None,
        cpu_equi_step: bool= False,
        cpu_equi_job_config: Dict= None,
        ) -> Dict:
    """this function prepare files for a submission ready MD task of HPCs.
    The task is to perform a production run of molecular dynamics simulation with the
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
    Returns:
        a dictionary in the structure of
            {
            "structure_files" : List[str],
            "job_list" : List[ClusterJob],
            }
    """
    # 1. make building blocks
    param_method, (
        min_step, heat_step, equi_step_1, equi_step_2, prod_step
    ) = _process_equi_md_sampling_arguments(
        stru = stru,
        param_method = param_method,
        prod_time = prod_time,
        prod_temperature = prod_temperature,
        prod_constrain = prod_constrain,
        record_period = record_period,
        cluster_job_config = cluster_job_config,
        cpu_equi_step = cpu_equi_step,
        cpu_equi_job_config = cpu_equi_job_config,
    )

    # 2. deploy files
    result = deployable_md_simulation(
        stru, param_method,
        steps=[min_step, heat_step, equi_step_1, equi_step_2, prod_step],
        parallel_runs=parallel_runs,
        work_dir=work_dir,
        )

    return result

def _process_equi_md_sampling_arguments(
        stru: Structure,
        param_method: MolDynParameterizer,
        # config for steps
        prod_time: float= 50.0, # ns
        prod_temperature: float = 300.0, #K
        prod_constrain: List[stru_cons.StructureConstraint]= None,
        record_period: float= 0.5, # ns
        cluster_job_config: Dict= None,
        cpu_equi_step: bool= False,
        cpu_equi_job_config: Dict= None,
        dont_freeze_bb_in_min: bool = False,
    ):
    """process the arguments of equi_md_sampling and deployable_equi_md_sampling
    into MolDynStep()s and MolDynParameterizer()"""
    if prod_constrain is None:
        prod_constrain = []
    for cons in prod_constrain:
        if cons.topology is not stru:
            # TODO probably can also support a eq-like checker in Structure()
            # that checks for the same topology and indexing.
            _LOGGER.error("inconsistency between constraint and stru. "
                          f"constraint topology: {cons.topology} "
                          f"stru: {stru} ")
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

    freeze_backbone = stru_cons.create_backbone_freeze(stru)
    # allow lifting the constriant in min so that backbone are relaxed upon extreme mutations
    if dont_freeze_bb_in_min:
        min_constrain = prod_constrain
    else:
        min_constrain = [freeze_backbone] + prod_constrain

    min_step  = parent_interface.build_md_step(
        name="min_micro",
        minimize=True,
        length=20000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=min_constrain)

    heat_step = parent_interface.build_md_step(
        name="heat_nvt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=[(0, 0), (0.05*0.9, prod_temperature), (-1, prod_temperature)],
        pressure_scaling="none",
        constrain=[freeze_backbone] + prod_constrain)

    equi_step_1 = parent_interface.build_md_step(
        name="equi_npt",
        length=prod_time * 0.01,
        cluster_job_config=equi_job_config,
        core_type=equi_core,
        temperature=prod_temperature,
        constrain=[freeze_backbone] + prod_constrain)

    equi_step_2 = parent_interface.build_md_step(
        name="equi_npt_free_bb",
        length=prod_time * 0.01,
        cluster_job_config=equi_job_config,
        core_type=equi_core,
        temperature=prod_temperature,
        constrain=prod_constrain)

    prod_step = parent_interface.build_md_step(
        name="prod_npt",
        length=prod_time,
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        restart=True,
        if_report=True,
        temperature=prod_temperature,
        record_period=record_period,
        constrain=prod_constrain)

    return param_method, [min_step, heat_step, equi_step_1, equi_step_2, prod_step]

# == general building blocks ==
def md_simulation(stru: Structure,
                  param_method: MolDynParameterizer,
                  steps: List[MolDynStep],
                  params_in: MolDynParameter = None,
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
    #   - params
    if not isinstance(param_method, MolDynParameterizer) and not isinstance(params_in, MolDynParameter):
        _LOGGER.error(f"Please either provide param_method (current: {param_method}) or params_in (current: {params_in})")
        raise ValueError
    if isinstance(param_method, MolDynParameterizer) and isinstance(params_in, MolDynParameter):
        _LOGGER.warning("Both param_method and params_in are provided. Only params_in is used.")
    if isinstance(params_in, MolDynParameter):
        params_engine = params_in.engine
    else:
        params_engine = param_method.engine
    #   - MD engine consistency
    for i, step in enumerate(steps):
        if step.engine != params_engine:
            _LOGGER.error(
                f"The engine of step #{i} ({step.engine}) does not match the parameterizer/parameter ({params_engine})!")
            raise InconsistentMDEngine
    #   - Parallel method suppoort
    if parallel_method and (parallel_method not in supported_parallel_method):
        _LOGGER.error(
            f"Parallel method: {parallel_method} is not in the supported list: {supported_parallel_method}")
        raise ValueError

    # II. make work dir
    fs.safe_mkdir(work_dir)

    # III. parameterize
    if isinstance(params_in, MolDynParameter):
        params = params_in
    else:
        param_method = copy.deepcopy(param_method)
        param_method.parameterizer_temp_dir = work_dir
        params = param_method.run(stru)

    # IV. run MD steps
    ## parallelize
    if parallel_method == "cluster_job":
        results = _parallelize_md_steps_with_cluster_job(parallel_runs, work_dir, steps, params, job_check_period)
    ## sequential
    if parallel_method is None:
        results = _serial_md_steps(parallel_runs, work_dir, steps, params)

    return params, results

def deployable_md_simulation(
        stru: Structure,
        param_method: MolDynParameterizer,
        steps: List[MolDynStep],
        parallel_runs: int=1,
        work_dir: str="./MD",
    ) -> Dict[str, List]:
    """This science API deploy a Molecular Dynamics simulation task as submission
    ready files.

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
        work_dir:
            the directory that contains all the MD files input/intermediate/output

    Return:
        a dictionary in the structure below
            {
            "structure_files" : List[str],
            "job_list" : List[ClusterJob],
            }"""
    # I. san check
    #   - MD engine consistency
    for i, step in enumerate(steps):
        if step.engine != param_method.engine:
            _LOGGER.error(
                f"The engine of step #{i} ({step.engine}) does not match the parameterizer ({param_method.engine})!")
            raise InconsistentMDEngine

    # II. make work dir
    fs.safe_mkdir(work_dir)

    # III. parameterize
    param_method = copy.deepcopy(param_method)
    param_method.parameterizer_temp_dir = work_dir
    params = param_method.run(stru)

    # IV. generate MD files
    results = []
    for i in range(parallel_runs):
        job_list = []
        result_egg_ele = []
        # create job path
        sub_work_dir = fs.get_valid_temp_name(f"{work_dir}/rep_{i}")
        fs.safe_mkdir(sub_work_dir)
        output = None  # the output place holder; the output between steps are very different for different packages so it will prob also becomes a class

        for step in steps:
            # 1. make job list
            step.work_dir = sub_work_dir # NOTE(QZ): this could case a bug. Consider deepcopy it
            if output: # steps after will use output from the previous one
                job, output = step.make_job(output, path_rel_to=sub_work_dir)
            else: # the 1st step
                job, output = step.make_job(params, path_rel_to=sub_work_dir)
            job_list.append(job)
            # 2. make output (we need to translate all since error checking and cleaning is needed)
            result_egg_ele.append((step, output))

        job_list = type(step).try_merge_jobs(job_list)

        parallel_result = {
            "structure_files" : params.file_list,
            "job_list" : job_list,
        }
        results.append(parallel_result)

    return results

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
        sub_work_dir = fs.get_valid_temp_name(f"{work_dir}/rep_{i}")
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
        sub_work_dir = fs.get_valid_temp_name(f"{work_dir}/rep_{i:06d}")
        fs.safe_mkdir(sub_work_dir)
        output = None
        result_ele = []
        for step in steps:
            step.work_dir = sub_work_dir
            if not output:
                output = step.run(params)
            else:
                output = step.run(output)
            result_ele.append(output)

        results.append(result_ele)

    return results

# ---------------------------------------------------------------------------
# Umbrella sampling Science API
# ---------------------------------------------------------------------------

def umbrella_sampling(
    stru: Structure,
    cv,  # CollectiveVariable – typed as Any to avoid circular import at module level
    window_targets,  # CVTargets | List[float]
    param_method: MolDynParameterizer = None,
    params_in: MolDynParameter = None,
    parallel_runs: int = 1,
    parallel_method: str = "cluster_job",
    work_dir: str = "./umbrella_MD",
    prod_time: float = 20.0,
    prod_temperature: float = 300.0,
    record_period: float = 0.2,
    cluster_job_config: Optional[Dict] = None,
    extra_constrain: Optional[List] = None,
    steps_builder: Optional[Callable] = None,
    job_check_period: int = 210,
) -> UmbrellaSamplingResult:
    """Run umbrella sampling across multiple CV windows.

    Each window is defined by a target value in *window_targets*.  For every
    window the function:

    1. Obtains per-window MD constraints from the CV (``AmberCV`` → via
       :py:meth:`~enzy_htp.structure.collective_variable.AmberCV.to_structure_constraint`).
    2. Optionally calls ``cv.serialize_for_engine(engine, window_dir, target)``
       to write engine-specific files and store the payload in
       :class:`UmbrellaWindowResult.metadata`.
    3. Builds MD steps via *steps_builder* (or a built-in equi+prod preset
       when *steps_builder* is ``None``).
    4. Invokes :func:`md_simulation` and collects results.

    Args:
        stru:
            Starting enzyme structure.
        cv:
            A :class:`~enzy_htp.structure.collective_variable.CollectiveVariable`
            instance (typically :class:`~enzy_htp.structure.collective_variable.AmberCV`).
        window_targets:
            An iterable of per-window RC target values (or a
            :class:`~enzy_htp.structure.collective_variable.CVTargets` object).
        param_method:
            MD parameterizer (e.g. ``AmberParameter``).  Mutually exclusive
            with *params_in* – providing *params_in* skips re-parameterization.
        params_in:
            Pre-computed :class:`MolDynParameter`.  When given, *param_method*
            is only used to determine the engine and build steps.
        parallel_runs:
            Number of parallel replicas per window.
        parallel_method:
            Parallelisation strategy passed to :func:`md_simulation`.
        work_dir:
            Root directory.  Per-window subdirectories are created as
            ``<work_dir>/window_<i>/``.
        prod_time:
            Production MD length (ns).  Used only when *steps_builder* is
            ``None``.
        prod_temperature:
            Production temperature (K).  Used only when *steps_builder* is
            ``None``.
        record_period:
            Trajectory recording period (ns).  Used only when *steps_builder*
            is ``None``.
        cluster_job_config:
            HPC job configuration dict forwarded to MD step builders.
        extra_constrain:
            Additional :class:`~enzy_htp.structure.structure_constraint.StructureConstraint`
            objects applied in every window (e.g. backbone freeze).
        steps_builder:
            Optional callable ``(constrain: List[StructureConstraint]) →
            List[MolDynStep]`` that fully controls step construction.  When
            supplied, *prod_time*, *prod_temperature*, *record_period*, and
            *cluster_job_config* are ignored (the builder should use them as
            needed via closure).
        job_check_period:
            Poll interval (s) for HPC job monitoring.

    Returns:
        A :class:`UmbrellaSamplingResult` containing per-window results and
        CV metadata.

    Raises:
        ValueError:    If neither *param_method* nor *params_in* is provided.
        NotImplementedError: If the CV type is not supported for automatic
                             constraint extraction and no *steps_builder* is
                             given.

    Example::

        from enzy_htp.structure import PDBParser
        from enzy_htp.structure.collective_variable import DistanceCV, AmberCV
        from enzy_htp.geometry import umbrella_sampling
        from enzy_htp._interface.amber_interface import AmberParameter

        stru  = PDBParser().get_structure("enzyme_amber.pdb")
        cv    = AmberCV(
                    DistanceCV("A.55.CA", "B.100.CA", name="reaction_coord"),
                    amber_params={"rk2": 100.0, "rk3": 100.0},
                )
        tgts  = cv.generate_window_targets(start=14.0, end=24.0, step=1.0)
        params = AmberParameter("enzyme.inpcrd", "enzyme.prmtop")
        result = umbrella_sampling(stru, cv, tgts, params_in=params,
                                   work_dir="./umbrella_MD",
                                   cluster_job_config={...})
    """
    # --- lazy import to avoid circular deps ---
    from enzy_htp.structure.collective_variable import AmberCV, CVTargets

    # --- sanity check ---
    if param_method is None and params_in is None:
        _LOGGER.error("umbrella_sampling: either param_method or params_in must be provided.")
        raise ValueError("Provide param_method or params_in.")

    if extra_constrain is None:
        extra_constrain = []

    # --- normalise window_targets ---
    if isinstance(window_targets, CVTargets):
        target_list: List[float] = list(window_targets.data)
    else:
        target_list = list(window_targets)

    # --- determine engine / interface ---
    if params_in is not None:
        engine = params_in.engine
        parent_interface = params_in.parent_interface if hasattr(params_in, "parent_interface") else None
    else:
        engine = param_method.engine
        parent_interface = param_method.parent_interface

    # Fall back to using param_method's interface when params_in doesn't expose it
    if parent_interface is None and param_method is not None:
        parent_interface = param_method.parent_interface

    # --- root work dir ---
    fs.safe_mkdir(work_dir)

    # --- collect results ---
    window_results: List[UmbrellaWindowResult] = []

    for i, target in enumerate(target_list):
        window_dir = fs.get_valid_temp_name(f"{work_dir}/window_{i:04d}")
        fs.safe_mkdir(window_dir)

        _LOGGER.info(
            f"umbrella_sampling: window {i}/{len(target_list) - 1} "
            f"target={target} {cv.unit} dir={window_dir}"
        )

        # ---- 1. Build per-window CV constraint ----
        cv_constraint = None
        if isinstance(cv, AmberCV):
            cv_constraint = cv.to_structure_constraint(
                stru, target, rs_filepath="{mdstep_dir}/cv.rs"
            )
        elif steps_builder is None:
            raise NotImplementedError(
                f"umbrella_sampling: automatic constraint extraction is only "
                f"supported for AmberCV, got '{type(cv).__name__}'.  "
                "Provide a steps_builder callable for other CV types."
            )

        # ---- 2. Assemble full constraint list ----
        window_constrain = (
            ([cv_constraint] if cv_constraint is not None else []) + extra_constrain
        )

        # ---- 3. Build MD steps ----
        if steps_builder is not None:
            steps = steps_builder(window_constrain)
        else:
            if parent_interface is None:
                raise ValueError(
                    "umbrella_sampling: cannot build default MD steps without "
                    "a parent_interface.  Provide a steps_builder or ensure "
                    "param_method exposes parent_interface."
                )
            steps = _build_default_umbrella_steps(
                parent_interface,
                stru,
                window_constrain,
                prod_time=prod_time,
                prod_temperature=prod_temperature,
                record_period=record_period,
                cluster_job_config=cluster_job_config,
            )

        # ---- 4. Run MD simulation ----
        _, md_result = md_simulation(
            stru=stru,
            param_method=param_method,
            steps=steps,
            params_in=params_in,
            parallel_runs=parallel_runs,
            parallel_method=parallel_method,
            work_dir=window_dir,
            job_check_period=job_check_period,
        )

        # ---- 5. Collect engine payload for metadata ----
        try:
            engine_payload = cv.serialize_for_engine(engine, window_dir, target)
        except (NotImplementedError, Exception) as exc:
            _LOGGER.debug(
                f"umbrella_sampling: cv.serialize_for_engine skipped for window {i}: {exc}"
            )
            engine_payload = {}

        window_results.append(UmbrellaWindowResult(
            window_index=i,
            window_target=target,
            md_results=md_result,
            metadata={
                "cv": cv.to_dict(),
                "engine_payload": engine_payload,
                "work_dir": window_dir,
            },
        ))

    return UmbrellaSamplingResult(
        windows=window_results,
        cv_dict=cv.to_dict(),
        window_targets=target_list,
        work_dir=work_dir,
    )


def _build_default_umbrella_steps(
    parent_interface,
    stru: Structure,
    constrain: List,
    prod_time: float = 20.0,
    prod_temperature: float = 300.0,
    record_period: float = 0.2,
    cluster_job_config: Optional[Dict] = None,
) -> List[MolDynStep]:
    """Build a default min → equi → prod MD step sequence for umbrella sampling.

    Args:
        parent_interface:    MD engine interface (e.g. ``interface.amber``).
        stru:                Structure for backbone-freeze generation.
        constrain:           Per-window constraints (CV + any extras).
        prod_time:           Production MD length (ns).
        prod_temperature:    Production temperature (K).
        record_period:       Trajectory recording period (ns).
        cluster_job_config:  HPC job configuration.

    Returns:
        ``[min_step, equi_step, prod_step]``
    """
    bb_freeze = stru_cons.create_backbone_freeze(stru)
    min_constrain = [bb_freeze] + constrain

    min_step = parent_interface.build_md_step(
        name="min_micro",
        minimize=True,
        length=20000,  # cycles
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=min_constrain,
    )

    equi_step = parent_interface.build_md_step(
        name="equi_npt",
        length=prod_time * 0.01,
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=prod_temperature,
        constrain=constrain,
    )

    prod_step = parent_interface.build_md_step(
        name="prod_npt",
        length=prod_time,
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        restart=True,
        if_report=True,
        temperature=prod_temperature,
        record_period=record_period,
        constrain=constrain,
    )

    return [min_step, equi_step, prod_step]


# == helper tools ==
def get_deployable_md_cli() -> str:
    """get the content of a CLI tool that manage deployed
    MD tasks in batch"""

