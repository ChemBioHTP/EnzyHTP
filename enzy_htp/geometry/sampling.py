"""Define functions for gemoetry sampling of Structure(). These functions sample gemoetries of
the target start structure on the energy surface defined by a certain energy function.

Science API:
    + md_simulation()
    + equi_md_sampling()
    + deployable_equi_md_sampling()
    + md_energy_injection()

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-7-30
"""
import copy
import math
import os
import sys
from typing import List, Dict, Tuple

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import InconsistentMDEngine
from enzy_htp.core import job_manager
from enzy_htp.core.general import save_obj, load_obj, save_func_to_main
from enzy_htp.structure import Structure, StructureEnsemble, StruSelection
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp.structure.structure_selection import select_stru
from enzy_htp._interface.handle_types import (
    MolDynStep,
    MolDynParameterizer,
    MolDynParameter,
    MolDynResult)


def _stitch_energy_injection_replica_trajectories(
        replica_state: Dict,
        parent_interface,
        work_dir: str,
    ) -> None:
    """Combine the production segment trajectories from one energy-injection replica.

    This helper collects the ordered production `.nc` files recorded for a single
    segmented local-heating replica and uses the engine interface to generate:

    1. an auto-imaged combined NetCDF trajectory

    The output filenames follow the default convention used by
    `md_energy_injection()`:

    - `{work_dir}/prod_npt_combined.nc`

    Args:
        replica_state:
            the serialized state dictionary for one segmented energy-injection replica.
            The state must contain `production_results` with `traj_file` entries and
            a `topology_file`.
        parent_interface:
            the engine interface used to combine segment trajectories. In the current
            implementation this is the Amber interface.
        work_dir:
            the per-replica working directory that will receive the combined output files.
    """
    combined_traj_path = os.path.join(work_dir, "prod_npt_combined.nc")
    replica_state["combined_traj_file"] = combined_traj_path

    if os.path.exists(combined_traj_path):
        replica_state["stitched_trajectory"] = True
        return

    traj_paths = []
    for result_record in replica_state["production_results"]:
        traj_path = result_record.get("traj_file")
        if traj_path is None:
            _LOGGER.warning(
                "Skipping stitched trajectory creation because a production segment has no traj_file."
            )
            replica_state["combined_traj_file"] = None
            replica_state["stitched_trajectory"] = False
            return
        traj_paths.append(traj_path)

    parent_interface.combine_traj_segments(
        traj_paths=traj_paths,
        topology_path=replica_state["topology_file"],
        out_path=combined_traj_path,
        autoimage=True,
    )
    replica_state["stitched_trajectory"] = True


def _build_energy_injection_md_steps(
        stru: Structure,
        param_method: MolDynParameterizer,
        segment_time: float,
        temperature: float,
        prod_constrain: List[stru_cons.StructureConstraint],
        record_period: float,
        cluster_job_config,
        cpu_equi_step: bool,
        cpu_equi_job_config,
        dont_freeze_bb_in_min: bool,
    ) -> Tuple[List[MolDynStep], MolDynStep]:
    """Build standard equilibration steps plus an ASCII-restart production template."""
    _, steps = _process_equi_md_sampling_arguments(
        stru=stru,
        param_method=param_method,
        prod_time=segment_time,
        prod_temperature=temperature,
        prod_constrain=prod_constrain,
        record_period=record_period,
        cluster_job_config=cluster_job_config,
        cpu_equi_step=cpu_equi_step,
        cpu_equi_job_config=cpu_equi_job_config,
        dont_freeze_bb_in_min=dont_freeze_bb_in_min,
    )
    eq_steps = list(steps[:-1])
    prod_step_template = steps[-1]
    prod_step_template.length = segment_time
    prod_step_template.restart = True
    prod_step_template.if_report = True
    prod_step_template.record_period = record_period
    if hasattr(prod_step_template, "ascii_rst"):
        prod_step_template.ascii_rst = True
    return eq_steps, prod_step_template


def _perturb_segment_restart(
        parent_interface,
        driven_selection: StruSelection,
        result_record: Dict,
        work_dir: str,
        segment_name: str,
        segment_idx: int,
        injection_mode: str,
        driving_temperature: float,
        remove_drift: bool,
    ) -> Tuple[str, Dict]:
    """Apply one inter-segment velocity perturbation and return the next restart."""
    modified_restart = os.path.join(work_dir, f"{segment_name}.modified.rst")
    segment_metadata = parent_interface.perturb_restart_velocities(
        selection=driven_selection,
        restart_in=result_record["last_frame_file"],
        restart_out=modified_restart,
        target_temperature=driving_temperature,
        mode=injection_mode,
        remove_drift=remove_drift,
    )
    segment_metadata["segment_index"] = segment_idx
    return modified_restart, segment_metadata


def _energy_injection_replica_state(
        checkpoint_path: str,
        params: MolDynParameter,
        segment_count: int,
        segment_time: float,
        prod_time: float,
    ) -> Dict:
    """Resume a compatible checkpoint or initialize a new replica state."""
    if checkpoint_path is None or not os.path.exists(checkpoint_path):
        return {
            "stage": "equilibration",
            "next_equilibration_step_idx": 0,
            "next_segment_idx": 0,
            "equilibration_results": [],
            "production_results": [],
            "segment_metadata": [],
            "current_restart": params.input_coordinate_file,
            "topology_file": params.topology_file,
            "segment_count": segment_count,
            "segment_time": segment_time,
            "prod_time": prod_time,
            "combined_traj_file": None,
            "stitched_trajectory": False,
        }

    state = load_obj(checkpoint_path)
    if state.get("segment_count") != segment_count:
        _LOGGER.error("Cannot resume energy-injection checkpoint with a different segment count.")
        raise ValueError
    saved_segment_time = state.get("segment_time")
    saved_prod_time = state.get("prod_time")
    if (
            saved_segment_time is None
            or saved_prod_time is None
            or not math.isclose(saved_segment_time, segment_time)
            or not math.isclose(saved_prod_time, prod_time)
    ):
        _LOGGER.error("Cannot resume energy-injection checkpoint with different production timing.")
        raise ValueError
    state.setdefault("combined_traj_file", None)
    state.setdefault("stitched_trajectory", False)
    return state


def _run_energy_injection_replica(
        stru: Structure,
        params: MolDynParameter,
        eq_steps: List[MolDynStep],
        prod_step_template: MolDynStep,
        driven_selection: StruSelection,
        segment_count: int,
        injection_mode: str,
        driving_temperature: float,
        prod_time: float,
        segment_time: float,
        work_dir: str,
        parent_interface,
        checkpoint_path: str,
        remove_drift: bool,
        stitched_trajectory: bool,
        replica_idx: int,
    ) -> Dict:
    """Run one complete local energy-injection replica."""
    fs.safe_mkdir(work_dir)
    state = _energy_injection_replica_state(
        checkpoint_path=checkpoint_path,
        params=params,
        segment_count=segment_count,
        segment_time=segment_time,
        prod_time=prod_time,
    )
    if checkpoint_path is not None:
        save_obj(state, checkpoint_path)

    current_restart = state["current_restart"]
    if state["next_equilibration_step_idx"] < len(eq_steps):
        eq_work_dir = os.path.join(work_dir, "equilibration")
        eq_params = params
        if state["equilibration_results"]:
            eq_params = type(params)(
                current_restart,
                params.topology_file,
                getattr(params, "ncaa_chrgspin_mapper", {}),
            )
        eq_results = []
        output = eq_params
        fs.safe_mkdir(eq_work_dir)
        for step in eq_steps[state["next_equilibration_step_idx"]:]:
            step.work_dir = eq_work_dir
            output = step.run(output)
            eq_results.append(output)
        eq_records = [
            {
                "traj_file": getattr(result, "traj_file", None),
                "traj_log_file": getattr(result, "traj_log_file", None),
                "last_frame_file": getattr(result, "last_frame_file", None),
                "source": getattr(result, "source", "amber"),
            }
            for result in eq_results
        ]
        state["equilibration_results"].extend(eq_records)
        state["next_equilibration_step_idx"] = len(eq_steps)
        state["current_restart"] = eq_records[-1]["last_frame_file"]
        current_restart = state["current_restart"]
        state["stage"] = "production"
        if checkpoint_path is not None:
            save_obj(state, checkpoint_path)

    for segment_idx in range(state["next_segment_idx"], segment_count):
        segment_step = copy.copy(prod_step_template)
        segment_step.length = min(segment_time, prod_time - (segment_idx * segment_time))
        segment_step.name = f"{prod_step_template.name}_seg_{segment_idx:06d}"
        segment_params = type(params)(
            current_restart,
            params.topology_file,
            getattr(params, "ncaa_chrgspin_mapper", {}),
        )
        segment_step.work_dir = os.path.join(work_dir, segment_step.name)
        fs.safe_mkdir(segment_step.work_dir)
        result = segment_step.run(segment_params)
        result_record = {
            "traj_file": getattr(result, "traj_file", None),
            "traj_log_file": getattr(result, "traj_log_file", None),
            "last_frame_file": getattr(result, "last_frame_file", None),
            "source": getattr(result, "source", "amber"),
        }
        state["production_results"].append(result_record)
        state["next_segment_idx"] = segment_idx + 1
        state["current_restart"] = result_record["last_frame_file"]
        current_restart = state["current_restart"]

        if segment_idx != segment_count - 1:
            current_restart, segment_metadata = _perturb_segment_restart(
                parent_interface=parent_interface,
                driven_selection=driven_selection,
                result_record=result_record,
                work_dir=work_dir,
                segment_name=segment_step.name,
                segment_idx=segment_idx,
                injection_mode=injection_mode,
                driving_temperature=driving_temperature,
                remove_drift=remove_drift,
            )
            state["segment_metadata"].append(segment_metadata)
            state["current_restart"] = current_restart

        state["stage"] = "complete" if state["next_segment_idx"] == segment_count else "production"
        if checkpoint_path is not None:
            save_obj(state, checkpoint_path)

    state["checkpoint_path"] = checkpoint_path
    if stitched_trajectory and state["stage"] == "complete":
        try:
            _stitch_energy_injection_replica_trajectories(
                replica_state=state,
                parent_interface=parent_interface,
                work_dir=work_dir,
            )
        except Exception as exc:
            _LOGGER.warning(
                f"Failed to stitch production trajectories for energy-injection replica {replica_idx}: {exc}"
            )
            state["combined_traj_file"] = None
            state["stitched_trajectory"] = False

    if checkpoint_path is not None:
        save_obj(state, checkpoint_path)
    return state


def md_energy_injection(
        stru: Structure,
        param_method: MolDynParameterizer,
        engine: str = "amber",
        work_dir: str = "./MD_SEGMENTED_INJECTION",
        prod_time: float = 10.0,
        segment_time: float = 0.05,
        temperature: float = 300.0,
        driven_region: str = None,
        injection_mode: str = "maxwell_reassign",
        driving_temperature: float = None,
        parallel_runs: int = 1,
        cluster_job_config: Dict = None,
        record_period: float = None,
        prod_constrain: List[stru_cons.StructureConstraint] = None,
        cpu_equi_step: bool = False,
        cpu_equi_job_config: Dict = None,
        dont_freeze_bb_in_min: bool = False,
        remove_drift: bool = True,
        checkpoint_fname: str = "energy_injection_checkpoint.pickle",
        result_fname: str = "energy_injection_result.pickle",
        stitched_trajectory: bool = True,
        ) -> Dict:
    """This science API performs segmented local-heating MD by repeatedly perturbing
    restart velocities in a selected region between short production segments.

    The workflow is:

    min (micro) -> heat (NVT) -> equi (NPT) -> equi free backbone (NPT)
    -> prod segment -> local velocity perturbation -> prod segment -> ...

    This API does **not** implement continuous local thermostatting during integration.
    Instead, it orchestrates standard Amber MD segments and mutates the restart
    velocities between segments in a user-selected driven region. This makes it useful
    for heuristic local-driving studies such as exploratory energy-transport analysis,
    qualitative response mapping, or repeated local perturbation tests.

    Supported perturbation modes are:

    - `maxwell_reassign`: redraw driven-region velocities from a Maxwell-Boltzmann
      distribution at `driving_temperature`
    - `velocity_scale`: rescale the existing driven-region velocities to match
      `driving_temperature`

    Args:
        stru:
            the starting structure
        param_method:
            the Parameterizer() used for parameterization. This determines the engine.
        engine:
            the molecular dynamics engine. Only `"amber"` is supported in the current
            implementation.
        work_dir:
            the directory that contains all the segmented-MD input/intermediate/output files
        prod_time:
            the total production simulation time across all segments (unit: ns)
        segment_time:
            the simulation time for each production segment between reinjection events
            (unit: ns)
        temperature:
            the base MD temperature used during equilibration and production
        driven_region:
            PyMOL-style selection string for the region that receives repeated local
            perturbations. Hydrogen atoms are always excluded.
        injection_mode:
            the local perturbation mode. Supported values are
            `"maxwell_reassign"` and `"velocity_scale"`.
        driving_temperature:
            the target temperature used in the driven region during the inter-segment
            perturbation step
        parallel_runs:
            the number of desired replicas
        cluster_job_config:
            the MD job config used to build the Amber MD steps
        record_period:
            the simulation time period for recording trajectory frames in each segment
            (unit: ns). If `None`, a default fraction of `segment_time` is used.
        prod_constrain:
            constraints applied in the production segments
        cpu_equi_step:
            whether to use CPU resources for equilibration steps
        cpu_equi_job_config:
            the job config for CPU equilibration if `cpu_equi_step=True`
        dont_freeze_bb_in_min:
            whether to remove the default backbone freeze during minimization
        remove_drift:
            whether to remove center-of-mass drift from the perturbed driven-region
            velocities after each inter-segment perturbation
        checkpoint_fname:
            filename used for the per-replica checkpoint state
        result_fname:
            filename used for the serialized workflow result saved under `work_dir`
        stitched_trajectory:
            whether to combine the production segment trajectories from each completed
            replica into an auto-imaged NetCDF trajectory (`prod_npt_combined.nc`)

    Returns:
        a dictionary in the structure of
            {
            "parameter" : MolDynParameter,
            "replicas" : List[Dict],
            "driven_region" : StruSelection,
            "segment_time" : float,
            "prod_time" : float,
            "injection_mode" : str,
            "driving_temperature" : float,
            }

        where each element of `replicas` contains the per-replica equilibration results,
        production results, segment metadata, checkpoint path, and when stitching is
        enabled/available:

            {
            ...
            "combined_traj_file" : str or None,
            "stitched_trajectory" : bool,
            }
    """
    if driven_region is None or driving_temperature is None:
        _LOGGER.error("driven_region and driving_temperature are required for md_energy_injection().")
        raise ValueError
    if injection_mode not in {"maxwell_reassign", "velocity_scale"}:
        _LOGGER.error(f"Unsupported injection_mode: {injection_mode}")
        raise ValueError
    if engine.lower() != "amber":
        _LOGGER.error(f"md_energy_injection() only supports Amber. Got engine={engine}")
        raise ValueError
    if param_method.engine.lower() != "amber":
        _LOGGER.error(f"Parameterizer engine {param_method.engine} is inconsistent with Amber energy injection.")
        raise InconsistentMDEngine

    segment_count = int(math.ceil(prod_time / segment_time))
    work_dir = os.path.abspath(work_dir)
    driven_atoms = [atom for atom in select_stru(stru, driven_region).atoms if not atom.is_hydrogen()]
    if not driven_atoms:
        _LOGGER.error(f"driven_region selects no non-hydrogen atoms: {driven_region}")
        raise ValueError
    driven_selection = StruSelection(driven_atoms)

    effective_record_period = record_period
    if effective_record_period is None:
        effective_record_period = segment_time * 0.1

    fs.safe_mkdir(work_dir)
    param_method = copy.deepcopy(param_method)
    param_method.parameterizer_temp_dir = work_dir
    params = param_method.run(stru)

    eq_steps, prod_step_template = _build_energy_injection_md_steps(
        stru=stru,
        param_method=param_method,
        segment_time=segment_time,
        temperature=temperature,
        prod_constrain=prod_constrain,
        cluster_job_config=cluster_job_config,
        record_period=effective_record_period,
        cpu_equi_step=cpu_equi_step,
        cpu_equi_job_config=cpu_equi_job_config,
        dont_freeze_bb_in_min=dont_freeze_bb_in_min,
    )

    result = {
        "parameter": params,
        "driven_region": driven_selection,
        "replicas": [],
    }
    saved_result = copy.copy(result)
    saved_result["replicas"] = []
    for replica_idx in range(parallel_runs):
        replica_dir = work_dir if parallel_runs == 1 else os.path.join(work_dir, f"rep_{replica_idx:06d}")
        checkpoint_path = os.path.join(replica_dir, checkpoint_fname)
        replica_state = _run_energy_injection_replica(
            stru=stru,
            params=params,
            eq_steps=eq_steps,
            prod_step_template=prod_step_template,
            driven_selection=driven_selection,
            segment_count=segment_count,
            injection_mode=injection_mode,
            driving_temperature=driving_temperature,
            prod_time=prod_time,
            segment_time=segment_time,
            work_dir=replica_dir,
            parent_interface=param_method.parent_interface,
            checkpoint_path=checkpoint_path,
            remove_drift=remove_drift,
            stitched_trajectory=stitched_trajectory,
            replica_idx=replica_idx,
        )
        topology_file = replica_state["topology_file"]
        result["replicas"].append({
            "equilibration": [
                param_method.parent_interface.deserialize_md_result(record, topology_file)
                for record in replica_state["equilibration_results"]
            ],
            "production": [
                param_method.parent_interface.deserialize_md_result(record, topology_file)
                for record in replica_state["production_results"]
            ],
            "injections": replica_state["segment_metadata"],
            "combined_traj_file": replica_state.get("combined_traj_file"),
            "checkpoint_path": checkpoint_path,
        })
        saved_result["replicas"].append({
            "equilibration": replica_state["equilibration_results"],
            "production": replica_state["production_results"],
            "injections": replica_state["segment_metadata"],
            "combined_traj_file": replica_state.get("combined_traj_file"),
            "checkpoint_path": checkpoint_path,
        })

    save_obj(saved_result, os.path.join(work_dir, result_fname))
    return result


def _energy_injection_deployable_main(
        sys_argvs,
        md_energy_injection_kwargs: dict,
    ):
    """Entry point for one submitted energy-injection replica workflow."""
    import os
    import enzy_htp.core.file_system as fs
    from enzy_htp.geometry import md_energy_injection

    result_path = None
    for idx, arg in enumerate(sys_argvs):
        if arg == "-o":
            result_path = sys_argvs[idx + 1]

    md_energy_injection_kwargs.pop("parallel_method", None)
    md_energy_injection_kwargs["parallel_runs"] = 1
    md_energy_injection(**md_energy_injection_kwargs)
    if result_path is not None:
        expected_result_path = os.path.abspath(os.path.join(
            md_energy_injection_kwargs["work_dir"],
            md_energy_injection_kwargs.get("result_fname", "energy_injection_result.pickle"),
        ))
        if os.path.abspath(result_path) != expected_result_path:
            fs.safe_cp(expected_result_path, result_path)


def deployable_md_energy_injection(
        stru: Structure,
        param_method: MolDynParameterizer,
        engine: str = "amber",
        work_dir: str = "./MD_SEGMENTED_INJECTION",
        cluster_job_config = None,
        main_fname: str = "energy_injection_main.py",
        kwargs_fname: str = "energy_injection_kwargs.pickle",
        result_fname: str = "energy_injection_result.pickle",
        sub_script_fname: str = "submit_energy_injection.sh",
        **kwargs,
    ) -> Dict[str, List]:
    """This function prepares files for a submission-ready segmented energy-injection
    MD task on HPCs.

    The generated task submits one GPU workflow job per replica. Each submitted job
    runs conventional equilibration, all production segments, and all inter-segment
    restart perturbations inside the same allocation.

    Args:
        stru:
            the starting structure
        param_method:
            the Parameterizer() used for parameterization. This determines the engine.
        engine:
            the molecular dynamics engine. Only `"amber"` is supported in the current
            implementation.
        work_dir:
            the directory that will contain generated replica workflow scripts, kwargs files,
            submission script, and eventual workflow outputs
        cluster_job_config:
            the cluster job configuration for the submitted segmented-MD replica worker jobs
        main_fname:
            filename for the generated Python main script in each replica directory
        kwargs_fname:
            filename for the serialized kwargs passed to the workflow main script
        result_fname:
            filename for the serialized workflow result
        sub_script_fname:
            filename for the generated cluster submission script
        **kwargs:
            additional keyword arguments forwarded to `md_energy_injection()`

    Returns:
        a dictionary in the structure of
            {
            "structure_files" : List[str],
            "job_list" : List[ClusterJob],
            "result_files" : List[str],
            "kwargs_files" : List[str],
            "main_files" : List[str],
            }
    """
    if engine.lower() != "amber":
        _LOGGER.error(f"Only Amber is supported in deployable_md_energy_injection(). Got {engine}")
        raise ValueError
    if param_method.engine.lower() != engine.lower():
        _LOGGER.error(
            f"The engine of the parameterizer ({param_method.engine}) does not match "
            f"the workflow engine ({engine})."
        )
        raise InconsistentMDEngine
    if cluster_job_config is None:
        _LOGGER.error("cluster_job_config is required for deployable_md_energy_injection().")
        raise ValueError

    work_dir = os.path.abspath(work_dir)
    if isinstance(cluster_job_config, job_manager.ClusterJobConfig):
        cluster_job_config = copy.deepcopy(cluster_job_config)
    elif isinstance(cluster_job_config, dict):
        cluster_job_config = job_manager.ClusterJobConfig.from_dict(cluster_job_config)
    else:
        _LOGGER.error(f"Unsupported cluster_job_config type: {type(cluster_job_config)}")
        raise TypeError
    if not cluster_job_config.has_cluster():
        _LOGGER.error("cluster_job_config must include a cluster for deployable_md_energy_injection().")
        raise ValueError

    fs.safe_mkdir(work_dir)
    cluster = cluster_job_config.cluster
    parallel_runs = kwargs.get("parallel_runs", 1)
    if parallel_runs < 1:
        _LOGGER.error(f"parallel_runs must be at least 1. Got {parallel_runs}")
        raise ValueError

    jobs = []
    result_paths = []
    kwargs_paths = []
    main_paths = []
    env_settings = {
        "head": getattr(cluster, "AMBER_ENV", {}).get("GPU", ""),
        "tail": "",
    }

    for replica_idx in range(parallel_runs):
        replica_dir = os.path.abspath(os.path.join(work_dir, f"rep_{replica_idx:06d}"))
        fs.safe_mkdir(replica_dir)
        main_path = os.path.abspath(os.path.join(replica_dir, main_fname))
        kwargs_path = os.path.abspath(os.path.join(replica_dir, kwargs_fname))
        result_path = os.path.abspath(os.path.join(replica_dir, result_fname))
        sub_script_path = os.path.abspath(os.path.join(replica_dir, sub_script_fname))

        replica_kwargs = dict(kwargs)
        replica_kwargs.update({
            "stru": stru,
            "param_method": param_method,
            "engine": engine,
            "work_dir": replica_dir,
            "parallel_runs": 1,
            "cluster_job_config": cluster_job_config,
            "result_fname": result_fname,
        })
        save_obj({"md_energy_injection_kwargs": replica_kwargs}, kwargs_path)
        save_func_to_main(_energy_injection_deployable_main, kwargs_path, main_path)

        res_keywords = copy.deepcopy(cluster_job_config.res_keywords)
        res_keywords.setdefault("job_name", f"md_energy_injection_{replica_idx:06d}")
        command = f"{os.path.abspath(sys.executable)} -u {main_path} -o {result_path} > {main_path}.out 2>&1"
        job = job_manager.ClusterJob.config_job(
            commands=command,
            cluster=cluster,
            env_settings=env_settings,
            res_keywords=res_keywords,
            sub_dir=replica_dir,
            sub_script_path=sub_script_path,
        )
        job.mimo = {
            "replica_idx": replica_idx,
            "replica_dir": replica_dir,
            "result_path": result_path,
            "kwargs_path": kwargs_path,
            "main_path": main_path,
        }
        jobs.append(job)
        result_paths.append(result_path)
        kwargs_paths.append(kwargs_path)
        main_paths.append(main_path)

    result = {
        "structure_files": [],
        "job_list": jobs,
        "result_files": result_paths,
        "kwargs_files": kwargs_paths,
        "main_files": main_paths,
    }
    if len(result_paths) == 1:
        result["result_file"] = result_paths[0]
        result["kwargs_file"] = kwargs_paths[0]
        result["main"] = main_paths[0]
    return result

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

# == helper tools ==
def get_deployable_md_cli() -> str:
    """get the content of a CLI tool that manage deployed
    MD tasks in batch"""
