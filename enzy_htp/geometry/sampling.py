"""Define functions for gemoetry sampling of Structure(). These functions sample gemoetries of
the target start structure on the energy surface defined by a certain energy function.

Science API:
    + md_simulation()
    + equi_md_sampling()
    + deployable_equi_md_sampling()

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-7-30
"""
import copy
import math
import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple
from collections.abc import Iterable

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import InconsistentMDEngine
from enzy_htp.core import job_manager
from enzy_htp.core.general import save_obj, load_obj, save_func_to_main
from enzy_htp.structure import Structure, StructureEnsemble, Atom, Residue, StruSelection
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp.structure.structure_selection import select_stru
from enzy_htp._interface.handle_types import (
    MolDynStep,
    MolDynParameterizer,
    MolDynParameter,
    MolDynResult)
_EXTREME_SEGMENT_WARNING_NS = 0.001
_EXTREME_DRIVING_TEMP_WARNING_K = 1000.0
_EXTREME_SEGMENT_COUNT_WARNING = 1000


def _normalize_region_to_stru_selection(
        stru: Structure,
        region,
        include_hydrogens: bool = False,
    ) -> StruSelection:
    """Normalize a user region specifier to a StruSelection."""
    if isinstance(region, StruSelection):
        atoms = list(region.atoms)
    elif isinstance(region, str):
        atoms = list(select_stru(stru, region).atoms)
    elif isinstance(region, Atom):
        atoms = [region]
    elif isinstance(region, Residue):
        atoms = list(region.atoms)
    elif isinstance(region, Iterable):
        atoms = []
        for item in region:
            if isinstance(item, Atom):
                atoms.append(item)
            elif isinstance(item, Residue):
                atoms.extend(item.atoms)
            else:
                _LOGGER.error(f"Unsupported region member type: {type(item)}")
                raise TypeError
    else:
        _LOGGER.error(f"Unsupported region type: {type(region)}")
        raise TypeError

    if not atoms:
        _LOGGER.error("Region selection is empty after normalization.")
        raise ValueError

    filtered_atoms = []
    for atom in atoms:
        if atom.root() is not stru:
            _LOGGER.error("Region atoms must belong to the supplied Structure.")
            raise ValueError
        if include_hydrogens or not atom.is_hydrogen():
            filtered_atoms.append(atom)

    if not filtered_atoms:
        _LOGGER.error("Region selection is empty after hydrogen filtering.")
        raise ValueError

    deduped_atoms = list(dict.fromkeys(filtered_atoms))
    return StruSelection(deduped_atoms)


def _serialize_md_result(result) -> Dict:
    """Convert an MD result-like object to a checkpoint-friendly record."""
    return {
        "traj_file": getattr(result, "traj_file", None),
        "traj_log_file": getattr(result, "traj_log_file", None),
        "last_frame_file": getattr(result, "last_frame_file", None),
        "source": getattr(result, "source", "amber"),
    }


def _stitch_energy_injection_replica_trajectories(
        replica_state: Dict,
        parent_interface,
        work_dir: str,
    ) -> None:
    """Combine the production segment trajectories from one energy-injection replica.

    This helper collects the ordered production `.nc` files recorded for a single
    segmented local-heating replica and uses the engine interface to generate:

    1. an auto-imaged combined NetCDF trajectory and
    2. a matching first-frame PDB file

    The output filenames follow the default convention used by
    `md_energy_injection()`:

    - `{work_dir}/prod_npt_combined.nc`
    - `{work_dir}/prod_npt_combined_frame1.pdb`

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
    combined_frame1_pdb_path = os.path.join(work_dir, "prod_npt_combined_frame1.pdb")
    replica_state["combined_traj_file"] = combined_traj_path
    replica_state["combined_frame1_pdb_file"] = combined_frame1_pdb_path

    if os.path.exists(combined_traj_path) and os.path.exists(combined_frame1_pdb_path):
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
            replica_state["combined_frame1_pdb_file"] = None
            replica_state["stitched_trajectory"] = False
            return
        traj_paths.append(traj_path)

    parent_interface.combine_traj_segments(
        traj_paths=traj_paths,
        topology_path=replica_state["topology_file"],
        out_path=combined_traj_path,
        frame1_pdb_path=combined_frame1_pdb_path,
        autoimage=True,
    )
    replica_state["stitched_trajectory"] = True


def _materialize_energy_injection_replica(replica_state: Dict, parent_interface) -> Dict:
    """Convert serialized per-step result records back to MolDynResult objects."""
    result = copy.deepcopy(replica_state)
    topology_file = replica_state["topology_file"]
    result["equilibration_results"] = [
        parent_interface.deserialize_md_result(record, topology_file)
        for record in replica_state["equilibration_results"]
    ]
    result["production_results"] = [
        parent_interface.deserialize_md_result(record, topology_file)
        for record in replica_state["production_results"]
    ]
    return result


def _checkpoint_energy_injection_state(state: Dict, checkpoint_path: str) -> None:
    """Persist energy-injection workflow state."""
    save_obj(state, checkpoint_path)


def _load_energy_injection_state(checkpoint_path: str) -> Dict:
    """Load persisted energy-injection workflow state."""
    if not os.path.exists(checkpoint_path):
        return None
    return load_obj(checkpoint_path)


def _coerce_cluster_job_config(cluster_job_config) -> job_manager.ClusterJobConfig:
    """Normalize cluster-job config to ClusterJobConfig."""
    if isinstance(cluster_job_config, job_manager.ClusterJobConfig):
        return copy.deepcopy(cluster_job_config)
    if isinstance(cluster_job_config, dict):
        return job_manager.ClusterJobConfig.from_dict(cluster_job_config)
    _LOGGER.error(f"Unsupported cluster_job_config type: {type(cluster_job_config)}")
    raise TypeError


def _get_runtime_python_env_settings(amber_env: str = "") -> Dict[str, str]:
    """Build child-job env settings from the current Python runtime.

    Prefer the active Amber runtime from the current process when present. Falling
    back to a cluster preset is only safe when the current environment does not
    already define Amber-related library paths.
    """
    current_python = os.path.abspath(sys.executable)
    current_prefix = os.path.dirname(os.path.dirname(current_python))
    current_path = os.environ.get("PATH", "")

    head_lines = [
        "unset PYTHONHOME",
        f"export PYTHONNOUSERSITE={os.environ.get('PYTHONNOUSERSITE', '1')}",
        f"export CONDA_PREFIX={current_prefix}",
    ]
    if current_path:
        head_lines.append(f"export PATH={current_path}")
    else:
        head_lines.append('export PATH="$CONDA_PREFIX/bin:$PATH"')
    current_pythonpath = os.environ.get("PYTHONPATH")
    if current_pythonpath:
        head_lines.append(f"export PYTHONPATH={current_pythonpath}")
    else:
        head_lines.append("unset PYTHONPATH")

    amber_runtime_vars = (
        "AMBERHOME",
        "CUDA_HOME",
        "LD_LIBRARY_PATH",
        "LIBRARY_PATH",
    )
    have_active_amber_runtime = False
    for env_var in amber_runtime_vars:
        env_value = os.environ.get(env_var)
        if env_value:
            head_lines.append(f"export {env_var}={env_value}")
            have_active_amber_runtime = True

    if amber_env and not have_active_amber_runtime:
        head_lines.append(amber_env)
    head_lines.append("hash -r")
    return {
        "head": os.linesep.join(head_lines),
        "tail": "",
    }


def _make_restart_parameter_from_existing(params: MolDynParameter, restart_path: str):
    """Reconstruct a parameter-like input from an existing restart path."""
    return type(params)(
        restart_path,
        params.topology_file,
        getattr(params, "ncaa_chrgspin_mapper", {}),
    )


def _run_segmented_local_replica(
        params: MolDynParameter,
        eq_steps: List[MolDynStep],
        prod_step_template: MolDynStep,
        driven_selection: StruSelection,
        injection_mode: str,
        driving_temperature: float,
        prod_time: float,
        segment_time: float,
        work_dir: str,
        parent_interface,
        checkpoint_path: str = None,
        remove_drift: bool = True,
    ) -> Dict:
    """Run one segmented local replica with inter-segment restart mutation."""
    segment_count = int(math.ceil(prod_time / segment_time))
    state = None
    if checkpoint_path is not None:
        state = _load_energy_injection_state(checkpoint_path)

    if state is None:
        state = {
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
        }
        if checkpoint_path is not None:
            _checkpoint_energy_injection_state(state, checkpoint_path)

    current_restart = state["current_restart"]

    for eq_idx in range(state["next_equilibration_step_idx"], len(eq_steps)):
        # Keep the original constraint/selection objects attached to the live topology.
        # Deep-copying step constraints can strand copied Atom objects without their
        # parent Structure, which breaks Amber index lookup.
        local_step = copy.copy(eq_steps[eq_idx])
        local_step.work_dir = work_dir
        if eq_idx == 0 and state["next_equilibration_step_idx"] == 0:
            current_input = params
        else:
            current_input = _make_restart_parameter_from_existing(params, current_restart)
        result = local_step.run(current_input)
        serialized = _serialize_md_result(result)
        state["equilibration_results"].append(serialized)
        state["next_equilibration_step_idx"] = eq_idx + 1
        state["current_restart"] = serialized["last_frame_file"]
        current_restart = serialized["last_frame_file"]
        if state["next_equilibration_step_idx"] == len(eq_steps):
            state["stage"] = "production"
        if checkpoint_path is not None:
            _checkpoint_energy_injection_state(state, checkpoint_path)

    for segment_idx in range(state["next_segment_idx"], segment_count):
        local_step = copy.copy(prod_step_template)
        local_step.work_dir = work_dir
        remaining_time = prod_time - (segment_idx * segment_time)
        local_step.length = min(segment_time, remaining_time)
        local_step.name = f"{prod_step_template.name}_seg_{segment_idx:06d}"

        segment_input = _make_restart_parameter_from_existing(params, current_restart)
        segment_result = local_step.run(segment_input)
        serialized = _serialize_md_result(segment_result)
        state["production_results"].append(serialized)
        state["next_segment_idx"] = segment_idx + 1
        state["current_restart"] = serialized["last_frame_file"]
        current_restart = serialized["last_frame_file"]

        if segment_idx != segment_count - 1:
            modified_restart = os.path.join(work_dir, f"{local_step.name}.modified.rst")
            segment_metadata = parent_interface.perturb_restart_velocities(
                selection=driven_selection,
                restart_in=serialized["last_frame_file"],
                restart_out=modified_restart,
                target_temperature=driving_temperature,
                mode=injection_mode,
                remove_drift=remove_drift,
            )
            segment_metadata["segment_index"] = segment_idx
            state["segment_metadata"].append(segment_metadata)
            state["current_restart"] = modified_restart
            current_restart = modified_restart

        state["stage"] = "complete" if state["next_segment_idx"] == segment_count else "production"
        if checkpoint_path is not None:
            _checkpoint_energy_injection_state(state, checkpoint_path)

    return state


def _finalize_energy_injection_replica_state(
        replica_state: Dict,
        checkpoint_path: str,
        stitched_trajectory: bool,
        parent_interface,
        work_dir: str,
        replica_idx: int = None,
    ) -> Dict:
    """Attach derived metadata to one serialized replica state."""
    replica_state.setdefault("combined_traj_file", None)
    replica_state.setdefault("combined_frame1_pdb_file", None)
    replica_state.setdefault("stitched_trajectory", False)
    if stitched_trajectory and replica_state["stage"] == "complete":
        try:
            _stitch_energy_injection_replica_trajectories(
                replica_state=replica_state,
                parent_interface=parent_interface,
                work_dir=work_dir,
            )
        except Exception as exc:
            rep_label = "" if replica_idx is None else f" {replica_idx}"
            _LOGGER.warning(
                f"Failed to stitch production trajectories for energy-injection replica{rep_label}: {exc}"
            )
            replica_state["combined_traj_file"] = None
            replica_state["combined_frame1_pdb_file"] = None
            replica_state["stitched_trajectory"] = False
    replica_state["checkpoint_path"] = checkpoint_path
    return replica_state


def _run_segmented_local_workflow(
        params: MolDynParameter,
        eq_steps: List[MolDynStep],
        prod_step_template: MolDynStep,
        driven_selection: StruSelection,
        injection_mode: str,
        driving_temperature: float,
        prod_time: float,
        segment_time: float,
        work_dir: str,
        parallel_runs: int,
        checkpoint_fname: str,
        remove_drift: bool,
        stitched_trajectory: bool,
        parent_interface,
    ) -> Tuple[List[Dict], List[Dict]]:
    """Run all replicas locally and return both live and serialized result records."""
    replicas = []
    serializable_replicas = []
    for replica_idx in range(parallel_runs):
        replica_dir = f"{work_dir}/rep_{replica_idx:06d}"
        fs.safe_mkdir(replica_dir)
        checkpoint_path = os.path.join(replica_dir, checkpoint_fname)
        replica_state = _run_segmented_local_replica(
            params=params,
            eq_steps=eq_steps,
            prod_step_template=prod_step_template,
            driven_selection=driven_selection,
            injection_mode=injection_mode,
            driving_temperature=driving_temperature,
            prod_time=prod_time,
            segment_time=segment_time,
            work_dir=replica_dir,
            parent_interface=parent_interface,
            checkpoint_path=checkpoint_path,
            remove_drift=remove_drift,
        )
        replica_state = _finalize_energy_injection_replica_state(
            replica_state=replica_state,
            checkpoint_path=checkpoint_path,
            stitched_trajectory=stitched_trajectory,
            parent_interface=parent_interface,
            work_dir=replica_dir,
            replica_idx=replica_idx,
        )
        serializable_replicas.append(replica_state)
        replicas.append(_materialize_energy_injection_replica(replica_state, parent_interface))

    return replicas, serializable_replicas


def _energy_injection_replica_child_main(
        sys_argvs,
        replica_child_job_kwargs: dict,
    ):
    """Child main script for one submitted segmented-MD replica worker job."""
    import os
    from enzy_htp.core.general import save_obj
    from enzy_htp.geometry.sampling import (
        _finalize_energy_injection_replica_state,
        _run_segmented_local_replica,
    )

    result_path = replica_child_job_kwargs["result_path"]
    for idx, arg in enumerate(sys_argvs):
        if arg == "-o":
            result_path = sys_argvs[idx + 1]

    replica_state = _run_segmented_local_replica(**replica_child_job_kwargs["replica_run_kwargs"])
    replica_state = _finalize_energy_injection_replica_state(
        replica_state=replica_state,
        checkpoint_path=replica_child_job_kwargs["checkpoint_path"],
        stitched_trajectory=replica_child_job_kwargs["stitched_trajectory"],
        parent_interface=replica_child_job_kwargs["parent_interface"],
        work_dir=replica_child_job_kwargs["work_dir"],
        replica_idx=replica_child_job_kwargs.get("replica_idx"),
    )
    save_obj(replica_state, result_path)


def _make_energy_injection_replica_job(
        cluster_job_config: job_manager.ClusterJobConfig,
        replica_idx: int,
        replica_dir: str,
        replica_run_kwargs: Dict,
        stitched_trajectory: bool,
        parent_interface,
        child_main_fname: str = "energy_injection_replica_child_main.py",
        kwargs_fname: str = "energy_injection_replica_child_main_kwargs.pickle",
        result_fname: str = "energy_injection_replica_result.pickle",
        sub_script_fname: str = "submit_energy_injection_replica.sh",
    ) -> Tuple[job_manager.ClusterJob, str]:
    """Create a GPU worker ClusterJob for one segmented energy-injection replica."""
    cluster = cluster_job_config.cluster
    res_keywords = {
        "core_type": "gpu",
        "nodes": "1",
        "node_cores": "1",
        "job_name": f"md_energy_inj_rep_{replica_idx:06d}",
        "mem_per_core": "8G",
        "walltime": "3-00:00:00",
    }
    if cluster_job_config.has_res_keywords():
        res_keywords.update(cluster_job_config.res_keywords)

    fs.safe_mkdir(replica_dir)
    child_main_path = os.path.abspath(os.path.join(replica_dir, child_main_fname))
    kwargs_path = os.path.abspath(os.path.join(replica_dir, kwargs_fname))
    result_path = os.path.abspath(os.path.join(replica_dir, result_fname))
    sub_script_path = os.path.abspath(os.path.join(replica_dir, sub_script_fname))
    checkpoint_path = replica_run_kwargs["checkpoint_path"]

    save_obj(
        {
            "replica_child_job_kwargs": {
                "replica_run_kwargs": replica_run_kwargs,
                "checkpoint_path": checkpoint_path,
                "stitched_trajectory": stitched_trajectory,
                "parent_interface": parent_interface,
                "work_dir": replica_dir,
                "result_path": result_path,
                "replica_idx": replica_idx,
            }
        },
        kwargs_path,
    )
    save_func_to_main(_energy_injection_replica_child_main, kwargs_path, child_main_path)

    amber_gpu_env = getattr(cluster, "AMBER_ENV", {}).get("GPU", "")
    env_settings = _get_runtime_python_env_settings(amber_env=amber_gpu_env)
    command = f"{os.path.abspath(sys.executable)} -u {child_main_path} -o {result_path} > {child_main_path}.out 2>&1"
    job = job_manager.ClusterJob.config_job(
        commands=command,
        cluster=cluster,
        env_settings=env_settings,
        res_keywords=res_keywords,
        sub_dir=os.path.abspath(replica_dir),
        sub_script_path=sub_script_path,
    )
    job.mimo = {
        "replica_result_path": result_path,
        "replica_kwargs_path": kwargs_path,
        "child_main_path": child_main_path,
        "replica_dir": replica_dir,
        "replica_idx": replica_idx,
    }
    return job, result_path


def _run_segmented_cluster_replicas(
        params: MolDynParameter,
        eq_steps: List[MolDynStep],
        prod_step_template: MolDynStep,
        driven_selection: StruSelection,
        injection_mode: str,
        driving_temperature: float,
        prod_time: float,
        segment_time: float,
        work_dir: str,
        parallel_runs: int,
        cluster_job_config: job_manager.ClusterJobConfig,
        checkpoint_fname: str,
        remove_drift: bool,
        stitched_trajectory: bool,
        parent_interface,
        job_check_period: int,
    ) -> Tuple[List[Dict], List[Dict]]:
    """Submit one GPU worker job per replica and collect serialized results."""
    jobs = []
    result_paths = []
    for replica_idx in range(parallel_runs):
        replica_dir = f"{work_dir}/rep_{replica_idx:06d}"
        checkpoint_path = os.path.join(replica_dir, checkpoint_fname)
        replica_run_kwargs = {
            "params": params,
            "eq_steps": eq_steps,
            "prod_step_template": prod_step_template,
            "driven_selection": driven_selection,
            "injection_mode": injection_mode,
            "driving_temperature": driving_temperature,
            "prod_time": prod_time,
            "segment_time": segment_time,
            "work_dir": replica_dir,
            "parent_interface": parent_interface,
            "checkpoint_path": checkpoint_path,
            "remove_drift": remove_drift,
        }
        job, result_path = _make_energy_injection_replica_job(
            cluster_job_config=cluster_job_config,
            replica_idx=replica_idx,
            replica_dir=replica_dir,
            replica_run_kwargs=replica_run_kwargs,
            stitched_trajectory=stitched_trajectory,
            parent_interface=parent_interface,
        )
        jobs.append(job)
        result_paths.append(result_path)

    failed_jobs = job_manager.ClusterJob.wait_to_array_end(jobs, period=job_check_period)
    if failed_jobs:
        failed_job_ids = [job.job_id for job in failed_jobs]
        _LOGGER.error(f"Energy-injection replica jobs failed: {failed_job_ids}")
        raise RuntimeError

    serializable_replicas = [load_obj(result_path) for result_path in result_paths]
    replicas = [
        _materialize_energy_injection_replica(replica_state, parent_interface)
        for replica_state in serializable_replicas
    ]
    return replicas, serializable_replicas


def md_energy_injection(
        stru: Structure,
        param_method: MolDynParameterizer,
        engine: str = "amber",
        work_dir: str = "./MD_SEGMENTED_INJECTION",
        prod_time: float = 10.0,
        segment_time: float = 0.05,
        temperature: float = 300.0,
        driven_region = None,
        injection_mode: str = "maxwell_reassign",
        driving_temperature: float = None,
        include_hydrogens: bool = False,
        parallel_runs: int = 1,
        parallel_method: str = None,
        cluster_job_config: Dict = None,
        record_period: float = None,
        return_observables: bool = False,
        prod_constrain: List[stru_cons.StructureConstraint] = None,
        cpu_equi_step: bool = False,
        cpu_equi_job_config: Dict = None,
        dont_freeze_bb_in_min: bool = False,
        remove_drift: bool = True,
        job_check_period: int = 210,
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
            the region that receives repeated local perturbations. Supported input forms
            include a `StruSelection`, a pymol-style selection string, a `Residue`,
            an `Atom`, or an iterable of residues/atoms.
        injection_mode:
            the local perturbation mode. Supported values are
            `"maxwell_reassign"` and `"velocity_scale"`.
        driving_temperature:
            the target temperature used in the driven region during the inter-segment
            perturbation step
        include_hydrogens:
            whether hydrogens are included when normalizing `driven_region`
        parallel_runs:
            the number of desired replicas
        parallel_method:
            the method to parallelize multiple runs. Only `None` (local immediate
            execution in the current Python process) is implemented right now.
        cluster_job_config:
            the MD job config used to build the Amber MD steps
        record_period:
            the simulation time period for recording trajectory frames in each segment
            (unit: ns). If `None`, a default fraction of `segment_time` is used.
        return_observables:
            reserved for future use. Currently only workflow metadata/results are returned.
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
        job_check_period:
            reserved for future use in asynchronous/cluster-driven execution
        checkpoint_fname:
            filename used for the per-replica checkpoint state
        result_fname:
            filename used for the serialized workflow result saved under `work_dir`
        stitched_trajectory:
            whether to combine the production segment trajectories from each completed
            replica into an auto-imaged NetCDF trajectory (`prod_npt_combined.nc`) and
            a matching first-frame PDB (`prod_npt_combined_frame1.pdb`)

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
            "remove_drift" : bool,
            "stitched_trajectory" : bool,
            "return_observables" : bool,
            }

        where each element of `replicas` contains the per-replica equilibration results,
        production results, segment metadata, checkpoint path, and when stitching is
        enabled/available:

            {
            ...
            "combined_traj_file" : str or None,
            "combined_frame1_pdb_file" : str or None,
            "stitched_trajectory" : bool,
            }
    """
    if engine.lower() != "amber":
        _LOGGER.error(f"Only Amber is supported in md_energy_injection(). Got {engine}")
        raise ValueError
    if param_method.engine.lower() != engine.lower():
        _LOGGER.error(
            f"The engine of the parameterizer ({param_method.engine}) does not match "
            f"the workflow engine ({engine})."
        )
        raise InconsistentMDEngine
    if injection_mode not in {"maxwell_reassign", "velocity_scale"}:
        _LOGGER.error(
            f"Unsupported injection_mode '{injection_mode}'. "
            "Supported modes are 'maxwell_reassign' and 'velocity_scale'."
        )
        raise ValueError
    supported_parallel_methods = {None, "cluster_job"}
    if parallel_method not in supported_parallel_methods:
        _LOGGER.error(
            f"Parallel method: {parallel_method} is not in the supported list: {sorted(supported_parallel_methods, key=str)}"
        )
        raise ValueError
    if driven_region is None:
        _LOGGER.error("driven_region is required for md_energy_injection().")
        raise ValueError
    if driving_temperature is None:
        _LOGGER.error("driving_temperature is required for md_energy_injection().")
        raise ValueError
    if driving_temperature <= 0:
        _LOGGER.error(f"driving_temperature must be positive. Got {driving_temperature}")
        raise ValueError
    if segment_time <= 0:
        _LOGGER.error(f"segment_time must be positive. Got {segment_time}")
        raise ValueError
    if prod_time < segment_time:
        _LOGGER.error(f"prod_time must be >= segment_time. Got prod_time={prod_time}, segment_time={segment_time}")
        raise ValueError
    if parallel_runs < 1:
        _LOGGER.error(f"parallel_runs must be at least 1. Got {parallel_runs}")
        raise ValueError
    segment_count = int(math.ceil(prod_time / segment_time))
    if segment_count < 1:
        _LOGGER.error("The segmented production workflow requires at least one production segment.")
        raise ValueError
    if segment_time <= _EXTREME_SEGMENT_WARNING_NS:
        _LOGGER.warning(
            f"Segment length {segment_time} ns is extremely short and may cause heavy restart overhead."
        )
    if driving_temperature >= _EXTREME_DRIVING_TEMP_WARNING_K:
        _LOGGER.warning(
            f"Driving temperature {driving_temperature} K is very high and may destabilize the simulation."
        )
    if segment_count >= _EXTREME_SEGMENT_COUNT_WARNING:
        _LOGGER.warning(
            f"Segment count {segment_count} is very large and may require substantial walltime."
        )
    if return_observables:
        _LOGGER.warning("return_observables is not implemented yet; returning workflow metadata only.")

    work_dir = os.path.abspath(work_dir)

    normalized_cluster_job_config = None
    if parallel_method == "cluster_job":
        if cluster_job_config is None:
            _LOGGER.error("cluster_job_config is required when parallel_method='cluster_job'.")
            raise ValueError
        normalized_cluster_job_config = _coerce_cluster_job_config(cluster_job_config)
        if not normalized_cluster_job_config.has_cluster():
            _LOGGER.error("cluster_job_config must include a cluster when parallel_method='cluster_job'.")
            raise ValueError

    driven_selection = _normalize_region_to_stru_selection(
        stru=stru,
        region=driven_region,
        include_hydrogens=include_hydrogens,
    )

    effective_record_period = record_period
    if effective_record_period is None:
        effective_record_period = segment_time * 0.1

    fs.safe_mkdir(work_dir)
    param_method = copy.deepcopy(param_method)
    param_method.parameterizer_temp_dir = work_dir
    params = param_method.run(stru)

    _, (min_step, heat_step, equi_step_1, equi_step_2, prod_step_template) = _process_equi_md_sampling_arguments(
        stru=stru,
        param_method=param_method,
        prod_time=segment_time,
        prod_temperature=temperature,
        prod_constrain=prod_constrain,
        record_period=effective_record_period,
        cluster_job_config=cluster_job_config or "default",
        cpu_equi_step=cpu_equi_step,
        cpu_equi_job_config=cpu_equi_job_config,
        dont_freeze_bb_in_min=dont_freeze_bb_in_min,
    )

    eq_steps = [min_step, heat_step, equi_step_1, equi_step_2]
    prod_step_template.length = segment_time
    prod_step_template.restart = True
    prod_step_template.if_report = True
    prod_step_template.record_period = effective_record_period

    if parallel_method == "cluster_job":
        replicas, serializable_replicas = _run_segmented_cluster_replicas(
            params=params,
            eq_steps=eq_steps,
            prod_step_template=prod_step_template,
            driven_selection=driven_selection,
            injection_mode=injection_mode,
            driving_temperature=driving_temperature,
            prod_time=prod_time,
            segment_time=segment_time,
            work_dir=work_dir,
            parallel_runs=parallel_runs,
            cluster_job_config=normalized_cluster_job_config,
            checkpoint_fname=checkpoint_fname,
            remove_drift=remove_drift,
            stitched_trajectory=stitched_trajectory,
            parent_interface=param_method.parent_interface,
            job_check_period=job_check_period,
        )
    else:
        replicas, serializable_replicas = _run_segmented_local_workflow(
            params=params,
            eq_steps=eq_steps,
            prod_step_template=prod_step_template,
            driven_selection=driven_selection,
            injection_mode=injection_mode,
            driving_temperature=driving_temperature,
            prod_time=prod_time,
            segment_time=segment_time,
            work_dir=work_dir,
            parallel_runs=parallel_runs,
            checkpoint_fname=checkpoint_fname,
            remove_drift=remove_drift,
            stitched_trajectory=stitched_trajectory,
            parent_interface=param_method.parent_interface,
        )

    result = {
        "parameter": params,
        "replicas": replicas,
        "driven_region": driven_selection,
        "segment_time": segment_time,
        "prod_time": prod_time,
        "injection_mode": injection_mode,
        "driving_temperature": driving_temperature,
        "remove_drift": remove_drift,
        "stitched_trajectory": stitched_trajectory,
        "return_observables": return_observables,
    }
    save_obj(
        {
            "parameter": params,
            "replicas": serializable_replicas,
            "driven_region": driven_selection,
            "segment_time": segment_time,
            "prod_time": prod_time,
            "injection_mode": injection_mode,
            "driving_temperature": driving_temperature,
            "remove_drift": remove_drift,
            "stitched_trajectory": stitched_trajectory,
            "return_observables": return_observables,
        },
        os.path.join(work_dir, result_fname),
    )
    return result


def _energy_injection_child_main(
        sys_argvs,
        md_energy_injection_kwargs: dict,
    ):
    """Child main script for deployable single-job energy-injection runs."""
    import os
    import enzy_htp.core.file_system as fs
    from enzy_htp.geometry import md_energy_injection

    result_path = None
    for idx, arg in enumerate(sys_argvs):
        if arg == "-o":
            result_path = sys_argvs[idx + 1]

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
        child_main_fname: str = "energy_injection_child_main.py",
        kwargs_fname: str = "energy_injection_child_main_kwargs.pickle",
        result_fname: str = "energy_injection_result.pickle",
        sub_script_fname: str = "submit_energy_injection.sh",
        **kwargs,
    ) -> Dict[str, List]:
    """This function prepares files for a submission-ready segmented energy-injection
    MD task on HPCs.

    The generated task follows a manager/worker pattern:

    - submit one top-level manager job without GPU resources
    - inside that manager job, submit one GPU worker job per replica
    - inside each worker job, run all production segments for that replica locally
    - perform all inter-segment restart perturbations within the running worker job

    This avoids submitting one cluster job per production segment while also keeping
    the manager job free of direct GPU MD execution.

    Args:
        stru:
            the starting structure
        param_method:
            the Parameterizer() used for parameterization. This determines the engine.
        engine:
            the molecular dynamics engine. Only `"amber"` is supported in the current
            implementation.
        work_dir:
            the directory that will contain the generated child main script, kwargs file,
            submission script, and eventual workflow outputs
        cluster_job_config:
            the cluster job configuration for the submitted segmented-MD replica worker jobs
        child_main_fname:
            filename for the generated child Python main script
        kwargs_fname:
            filename for the serialized kwargs passed to the child main script
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
            "result_file" : str,
            "kwargs_file" : str,
            "child_main" : str,
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
    cluster_job_config = _coerce_cluster_job_config(cluster_job_config)
    if not cluster_job_config.has_cluster():
        _LOGGER.error("cluster_job_config must include a cluster for deployable_md_energy_injection().")
        raise ValueError

    cluster = cluster_job_config.cluster
    res_keywords = {
        "core_type": "cpu",
        "nodes": "1",
        "node_cores": "1",
        "job_name": "md_energy_injection_manager",
        "mem_per_core": "4G",
        "walltime": "3-00:00:00",
    }
    if cluster_job_config.has_res_keywords():
        manager_overrides = {
            k: v for k, v in cluster_job_config.res_keywords.items()
            if k in {"partition", "account", "qos", "exclude_node_id", "constraint", "walltime"}
        }
        res_keywords.update(manager_overrides)

    fs.safe_mkdir(work_dir)
    child_main_path = os.path.abspath(f"{work_dir}/{child_main_fname}")
    kwargs_path = os.path.abspath(f"{work_dir}/{kwargs_fname}")
    result_path = os.path.abspath(f"{work_dir}/{result_fname}")
    sub_script_path = os.path.abspath(f"{work_dir}/{sub_script_fname}")

    md_energy_injection_kwargs = {
        "stru": stru,
        "param_method": param_method,
        "engine": engine,
        "work_dir": work_dir,
        "parallel_method": "cluster_job",
        "cluster_job_config": cluster_job_config,
        "result_fname": result_fname,
    }
    md_energy_injection_kwargs.update(kwargs)
    save_obj({"md_energy_injection_kwargs": md_energy_injection_kwargs}, kwargs_path)
    save_func_to_main(_energy_injection_child_main, kwargs_path, child_main_path)

    amber_cpu_env = getattr(cluster, "AMBER_ENV", {}).get("CPU", "")
    env_settings = _get_runtime_python_env_settings(amber_env=amber_cpu_env)
    command = f"{os.path.abspath(sys.executable)} -u {child_main_path} -o {result_path} > {child_main_path}.out 2>&1"
    job = job_manager.ClusterJob.config_job(
        commands=command,
        cluster=cluster,
        env_settings=env_settings,
        res_keywords=res_keywords,
        sub_dir=os.path.abspath(work_dir),
        sub_script_path=sub_script_path,
    )

    return {
        "structure_files": [],
        "job_list": [job],
        "result_file": result_path,
        "kwargs_file": kwargs_path,
        "child_main": child_main_path,
    }

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
    if cluster_job_config is None:
        cluster_job_config = "default"

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
