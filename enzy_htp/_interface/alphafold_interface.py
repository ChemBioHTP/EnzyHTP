"""Defines an AlphafoldInterface class that serves as a bridge for enzy_htp to utilize AlphaFold2 software.

Author: Gemini, Claude Code
Date: 2025-08-16
"""
from __future__ import annotations
import tempfile
import os
import re
import copy
import json
from typing import Union, List, Optional, Dict, Tuple, Any
from pathlib import Path
from dataclasses import dataclass

from enzy_htp import config as eh_config
from .base_interface import BaseInterface
from .handle_types.modeling_engine import ModelingResultEgg
from enzy_htp.structure import Structure
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp._config.alphafold_config import AlphafoldConfig
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core.general import load_obj
from enzy_htp.chemical.sequence import create_fasta_from_sequences,parse_fasta_file


@dataclass
class AlphaFold2ResultEgg(ModelingResultEgg):
    """Result egg for AlphaFold2 predictions.
    
    Contains information about expected output files from AlphaFold2 predictions.
    Each result egg should only be associated with one job.
    """
    fasta_path: str
    """Path to the input FASTA file."""
    
    output_dir: Path
    """Directory where prediction results will be stored."""
    
    sequence_ids: List[str]
    """List of sequence identifiers in this job."""
    
    job: ClusterJob
    """Single cluster job for this result egg."""
    
    interface: 'AlphafoldInterface'
    """Reference to the parent interface for accessing helper methods."""


class AlphafoldInterface(BaseInterface):
    """Class that provides a direct interface for enzy_htp to utilize AlphaFold2 software."""

    def __init__(self, parent, config: AlphafoldConfig = None) -> None:
        """Simplistic constructor that optionally takes an AlphafoldConfig object as its only argument.
        Calls parent class."""
        super().__init__(parent, config, AlphafoldConfig)

    def af2_predict(
        self, 
        sequences: Union[List[str], List[List[str]]], 
        # Run type settings
        work_dir: Union[str, Path, None] = None,
        # -- local run
        non_armer_core_type: str = "gpu",
        # -- cluster job run
        cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None,
        array_size: int = 0,
        job_check_period: int = 30,
        seq_per_job: int = 1,
        # AlphaFold2 settings
        # -- model related
        model_preset: Optional[str] = None,
        num_models: int = 5,
        num_recycles: int = 3,
        num_multimer_predictions_per_model: int = 5, # native AF2 only
        # -- relax related
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        # -- msa related
        db_preset: str = "full_dbs", # native AF2 only
        use_precomputed_msas: bool = False, # native AF2 only
        # -- template related
        use_templates: bool = False,
        max_template_date: Optional[str] = None, # native AF2 only
        additional_options: Optional[List[str]] = None,
        random_seed: Optional[int] = None,
        **kwargs
    ) -> Dict[Union[str, Tuple[str, ...]], Dict[str, Any]]:
        """Science API for AlphaFold2 structure prediction.

        Args:
            sequences: Sequences to predict. Supported forms:
                      - List[str]: Multiple monomers [seq1, seq2, ...]
                      - List[List[str]]: Multimers per target [[A_seq, B_seq], [A_seq, B_seq, C_seq], ...]
            work_dir: Output/work directory. If None, creates a temporary directory.
            non_armer_core_type: Core type for local execution ("gpu" or "cpu").
            cluster_job_config: Configuration for cluster job submission.
            array_size: Number of jobs to run simultaneously (for cluster submission)
            job_check_period: Time cycle for job state checking (seconds)
            seq_per_job: Number of sequences per job (for array execution)
            num_models: Number of models to generate
            num_recycles: Number of recycling iterations
            num_multimer_predictions_per_model: Multimer predictions per model (native AF2)
            num_relax: Number of top ranked structures to relax
            relax_max_iteration: Maximum relaxation iterations
            db_preset: Database preset (native AF2)
            use_precomputed_msas: Whether to use precomputed MSAs (native AF2)
                (NOTE: when native AF2 distribution is used, HHSearch still runs.
                Details: https://github.com/google-deepmind/alphafold/issues/469)
            use_templates: Whether to use templates
            max_template_date: Maximum template date
            model_preset: Model preset configuration
            additional_options: Additional command-line options
                (NOTE: when native af2 distribution is used. HHSearch will still be run. Details: https://github.com/google-deepmind/alphafold/issues/469)
            random_seed: Random seed for reproducibility
            
        Returns:
            Dict mapping each original input to comprehensive prediction results.
            Keys are either the monomer sequence string, or a tuple[str, ...]
            for multimer inputs. Each value contains:
                - "best_model": Structure object of the best ranked model
                - "best_model_plddt": List of pLDDT scores for best model  
                - "best_model_index": Index (1-based) of the best model
                - "model_1": Structure object of model 1
                - "model_1_plddt": List of pLDDT scores for model 1
                - ... (continues for all models)
        """
        temp_paths = []
        stru_parser = PDBParser()
        if not isinstance(sequences, list):
            sequences = [sequences]
        # san check on sequence
        self._validate_sequences(sequences)
        
        try:
            # Create output directory first
            if work_dir is None:
                scratch_dir = eh_config.system.SCRATCH_DIR
                work_dir = fs.get_valid_temp_name(f"{scratch_dir}/alphafold2")
                temp_paths.append(work_dir)
            fs.safe_mkdir(work_dir)

            sequences_mapper = self._create_sequence_mapper(sequences, work_dir)

            if cluster_job_config:
                # Run on cluster using array jobs
                result_eggs = self.make_job(
                    sequences=sequences_mapper,
                    out_dir=work_dir,
                    num_models=num_models,
                    num_recycles=num_recycles,
                    num_relax=num_relax,
                    relax_max_iteration=relax_max_iteration,
                    use_templates=use_templates,
                    max_template_date=max_template_date,
                    model_preset=model_preset,
                    db_preset=db_preset,
                    additional_options=additional_options,
                    use_precomputed_msas=use_precomputed_msas,
                    random_seed=random_seed,
                    cluster_job_config=cluster_job_config,
                    seq_per_job=seq_per_job,
                    num_multimer_predictions_per_model=num_multimer_predictions_per_model
                )

                # Submit and wait for array jobs
                all_jobs = [egg.job for egg in result_eggs]
                failed_jobs = ClusterJob.wait_to_array_end_plus(
                    all_jobs,
                    period=job_check_period,
                    array_size=array_size
                )
                
                if failed_jobs:
                    _LOGGER.error(f"{len(failed_jobs)} AlphaFold2 jobs failed")
                    raise RuntimeError(f"{len(failed_jobs)} AlphaFold2 jobs failed")
            else:
                # Run locally
                self.run(
                    sequences=sequences_mapper,
                    out_dir=work_dir,
                    num_models=num_models,
                    num_recycles=num_recycles,
                    num_relax=num_relax,
                    relax_max_iteration=relax_max_iteration,
                    use_templates=use_templates,
                    max_template_date=max_template_date,
                    model_preset=model_preset,
                    db_preset=db_preset,
                    additional_options=additional_options,
                    use_precomputed_msas=use_precomputed_msas,
                    random_seed=random_seed,
                    non_armer_core_type=non_armer_core_type,
                    num_multimer_predictions_per_model=num_multimer_predictions_per_model
                )

            # Parse results into comprehensive output format
            comprehensive_results = {}
            for seq_id in sequences_mapper.keys():
                original_seq = sequences_mapper[seq_id]
                if isinstance(original_seq, list):
                    original_seq = tuple(original_seq)
                seq_results = self._parse_comprehensive_results(work_dir, seq_id, stru_parser)
                comprehensive_results[original_seq] = seq_results

            return comprehensive_results

        finally:
            # Clean up temporary files
            fs.clean_temp_file_n_dir(temp_paths)

    def run(
        self,
        sequences: Dict[str, Union[List[str], List[List[str]]]],
        out_dir: Union[str, Path],
        num_models: int = 5,
        num_recycles: int = 3,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = "alphafold2_ptm",
        db_preset: str = "full_dbs",
        additional_options: Optional[List[str]] = None,
        use_precomputed_msas: bool = False,
        random_seed: Optional[int] = None,
        non_armer_core_type: str = "gpu",
        num_multimer_predictions_per_model: int = 5,
    ) -> None:
        """Execute AlphaFold2 prediction locally.
        
        Args:
            sequences: Input sequences to predict {seq_id: [seq, ...], ...}
            out_dir: Output directory for results
            num_models: Number of models to generate
            num_recycles: Number of recycling iterations
            num_relax: Number of structures to relax
            relax_max_iteration: Maximum relaxation iterations
            use_templates: Whether to use templates
            max_template_date: Maximum template date
            model_preset: Model preset configuration
            db_preset: Database preset
            additional_options: Additional command-line options
            use_precomputed_msas: Whether to use precomputed MSAs
            random_seed: Random seed for reproducibility
            non_armer_core_type: Type of computing core ('gpu' or 'cpu') for local execution
            
        Returns:
            None. Executes the prediction; results are written to disk and
            parsed later by `_parse_comprehensive_results`.
        """
        config = self.config()
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)

        # Create FASTA file(s) from sequence data
        if self.config().INSTALL_TYPE != "colabfold_container":
            # Native AlphaFold with multimers: create separate FASTA files
            fasta_paths = self._native_af2_format_fasta(sequences, out_dir)
            fasta_path = ",".join(fasta_paths)
        else:
            # Single FASTA file for ColabFold
            fasta_path = self._colabfold_format_fasta(sequences, out_dir)
        
        # build command
        if config.INSTALL_TYPE == "colabfold_container":
            cmd = self._build_colabfold_container_command(
                fasta_path = fasta_path, 
                out_dir = out_dir, 
                num_models = num_models, 
                num_recycles = num_recycles, 
                num_relax = num_relax, 
                relax_max_iteration = relax_max_iteration, 
                use_templates = use_templates, 
                model_preset = model_preset, 
                additional_options = additional_options, 
                core_type = non_armer_core_type
            )
        elif config.INSTALL_TYPE == "alphafold2_native_container":
            cmd = self._build_alphafold2_native_container_command(
                fasta_path = fasta_path, 
                out_dir = out_dir, 
                num_relax = num_relax, 
                use_templates = use_templates,
                model_preset = model_preset,
                additional_options = additional_options,
                max_template_date = max_template_date,
                db_preset = db_preset,
                use_precomputed_msas = use_precomputed_msas,
                random_seed = random_seed,
                core_type = non_armer_core_type,
                num_multimer_predictions_per_model = num_multimer_predictions_per_model,
            )
        elif config.INSTALL_TYPE == "alphafold2_native_python":
            cmd = self._build_alphafold2_native_python_command(
                fasta_path = fasta_path, 
                out_dir = out_dir, 
                num_relax = num_relax, 
                use_templates = use_templates,
                model_preset = model_preset,
                additional_options = additional_options,
                max_template_date = max_template_date,
                db_preset = db_preset,
                use_precomputed_msas = use_precomputed_msas,
                random_seed = random_seed,
                core_type = non_armer_core_type,
                num_multimer_predictions_per_model = num_multimer_predictions_per_model,
            )
        else:
            _LOGGER.error(f"Unsupported install type: {config.INSTALL_TYPE}")
            raise ValueError(f"Unsupported install type: {config.INSTALL_TYPE}")
        
        _LOGGER.info(f"Running AlphaFold2 command: {' '.join(cmd)}")
        
        try:
            # Use env_manager to run command
            self.env_manager_.run_command(
                exe=cmd[0],
                args=cmd[1:],
            )
            _LOGGER.info(f"AlphaFold2 completed successfully")
            
        except Exception as e:
            _LOGGER.error(f"AlphaFold2 execution failed: {e}")
            raise RuntimeError(f"AlphaFold2 execution failed: {e}")

    def make_job(
        self,
        sequences: Dict[str, Union[List[str], List[List[str]]]],
        out_dir: Union[str, Path],
        num_models: int = 5,
        num_recycles: int = 3,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = None,
        db_preset: str = "full_dbs",
        additional_options: Optional[List[str]] = None,
        use_precomputed_msas: bool = False,
        random_seed: Optional[int] = None,
        cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None,
        seq_per_job: int = 1,
        num_multimer_predictions_per_model: int = 5
    ) -> List[AlphaFold2ResultEgg]:
        """Create cluster jobs for running AlphaFold2 with array support.
        
        Args:
            sequences: Mapping of sequence IDs to per-target sequences.
                Values are either
                - List[str]: monomer (single sequence)
                - List[List[str]]: multimer (list of chain sequences)
            out_dir: Output directory for results
            num_models: Number of models to generate
            num_recycles: Number of recycling iterations
            num_relax: Number of structures to relax
            relax_max_iteration: Maximum relaxation iterations
            use_templates: Whether to use templates
            max_template_date: Maximum template date
            model_preset: Model preset configuration
            db_preset: Database preset
            additional_options: Additional command-line options
            use_precomputed_msas: Whether to use precomputed MSAs
            random_seed: Random seed for reproducibility
            cluster_job_config: Configuration for cluster job submission
            seq_per_job: Number of sequences per job
            
        Returns:
            List of AlphaFold2ResultEgg objects, one per job
        """
        afconfig = self.config()
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)
        seq_w_id = list(sequences.items())
        
        # Create individual jobs for array execution
        result_eggs = []
        job_data = []
        for i in range(0, len(seq_w_id), seq_per_job):
            job_sequences_mapper = dict(seq_w_id[i:i+seq_per_job])
            if afconfig.INSTALL_TYPE == "colabfold_container":
                job_fasta_path = self._colabfold_format_fasta(job_sequences_mapper, out_dir)
            else:
                job_fasta_path = self._native_af2_format_fasta(job_sequences_mapper, out_dir)
                job_fasta_path = ",".join(job_fasta_path)
            job_data.append((job_fasta_path, list(job_sequences_mapper.keys())))
        # Create cluster jobs for all job data
        result_eggs = []
        
        for job_index, (job_fasta_path, job_seq_ids) in enumerate(job_data):
            
            # Create ClusterJob
            if isinstance(cluster_job_config, dict):
                job_config = ClusterJobConfig.from_dict(cluster_job_config)
            else:
                job_config = copy.deepcopy(cluster_job_config)
            
            # Use ClusterJob.config_job to create the job properly
            if not job_config.has_cluster():
                raise ValueError("cluster_job_config must specify a cluster for job execution")
            core_type = job_config.core_type if job_config.core_type else "gpu"
            
            # Build command for this job            
            if afconfig.INSTALL_TYPE == "colabfold_container":
                cmd = self._build_colabfold_container_command(
                    fasta_path=job_fasta_path,
                    out_dir=out_dir,
                    num_models=num_models,
                    num_recycles=num_recycles,
                    num_relax=num_relax,
                    relax_max_iteration=relax_max_iteration,
                    use_templates=use_templates,
                    model_preset=model_preset,
                    additional_options=additional_options,
                    core_type=core_type
                )
            elif afconfig.INSTALL_TYPE == "alphafold2_native_container":
                cmd = self._build_alphafold2_native_container_command(
                    fasta_path=job_fasta_path, 
                    out_dir=out_dir, 
                    num_relax=num_relax, 
                    use_templates=use_templates,
                    model_preset=model_preset, 
                    additional_options=additional_options,
                    max_template_date=max_template_date,
                    db_preset=db_preset,
                    use_precomputed_msas=use_precomputed_msas,
                    random_seed=random_seed,
                    core_type=core_type,
                    num_multimer_predictions_per_model=num_multimer_predictions_per_model
                )
            elif afconfig.INSTALL_TYPE == "alphafold2_native_python":
                cmd = self._build_alphafold2_native_python_command(
                    fasta_path=job_fasta_path, 
                    out_dir=out_dir, 
                    num_relax=num_relax, 
                    use_templates=use_templates,
                    model_preset=model_preset, 
                    additional_options=additional_options,
                    max_template_date=max_template_date,
                    db_preset=db_preset,
                    use_precomputed_msas=use_precomputed_msas,
                    random_seed=random_seed,
                    core_type=core_type,
                    num_multimer_predictions_per_model=num_multimer_predictions_per_model
                )
            else:
                raise ValueError(f"Unsupported install type: {afconfig.INSTALL_TYPE}")
            
            # Handle default res_keywords similar to amber_interface
            res_keywords_update = job_config.res_keywords if job_config.has_res_keywords() else {}
            default_res_keywords = afconfig.get_default_af2_cluster_job_res_keywords(core_type)
            job_config.res_keywords = default_res_keywords | res_keywords_update
            
            # Set up environment settings based on cluster and core type
            cluster = job_config.cluster
            if "container" in afconfig.INSTALL_TYPE:
                env_settings = cluster.CONTAINER_ENV[core_type.upper()]
            else:
                env_settings = cluster.AF2_ENV[core_type.upper()]
              
            # Create submission script path
            sub_script_path = out_dir / f"submit_alphafold_{job_index}.cmd"
            sub_script_path = fs.get_valid_temp_name(str(sub_script_path))
            
            job = ClusterJob.config_job(
                commands=' '.join(cmd),
                cluster=cluster,
                env_settings=env_settings,
                res_keywords=job_config.res_keywords,
                sub_dir="./",  # Use current dir for submission since paths are absolute
                sub_script_path=sub_script_path,
            )
            
            # Create result egg for this single job
            result_egg = AlphaFold2ResultEgg(
                fasta_path=job_fasta_path,
                output_dir=out_dir,
                sequence_ids=job_seq_ids,
                job=job,
                interface=self
            )
            
            result_eggs.append(result_egg)
        
        return result_eggs

    def _build_colabfold_container_command(
        self, 
        fasta_path: str, 
        out_dir: Path, 
        num_models: int,
        num_recycles: int,
        num_relax: int,
        relax_max_iteration: int,
        use_templates: bool,
        model_preset: str,
        additional_options: Optional[List[str]],
        core_type: str = "gpu"
    ) -> List[str]:
        """Build command for ColabFold container execution.
        
        Args:
            fasta_path: Path to input FASTA file.
            out_dir: Output directory (will be bind-mounted as /work).
            num_models: Number of models to generate.
            num_recycles: Number of recycling iterations.
            num_relax: Number of structures to relax with Amber.
            relax_max_iteration: Maximum iterations for relaxation.
            use_templates: Whether to enable template usage.
            model_preset: ColabFold/AF2 model preset to use.
            additional_options: Additional CLI options to append.
            core_type: Computing core type ('gpu' or 'cpu').
        """
        config = self.config_
        
        # Expand user paths
        container_path = os.path.expanduser(config.CONTAINER_PATH)
        
        cmd = [config.CONTAINER_TYPE, "run", "--cleanenv"]
        
        # Add GPU support if available and core type is GPU
        if core_type == "gpu" and config.CONTAINER_TYPE in ["docker", "apptainer", "singularity"]:
            cmd.append("--nv")
        
        # Add bind mounts
        for container_path_bind, host_path in config.CONTAINER_BIND_PATHS.items():
            expanded_host_path = os.path.expanduser(host_path)
            cmd.extend(["-B", f"{expanded_host_path}:{container_path_bind}"])
        
        # Bind work directory dynamically
        cmd.extend(["-B", f"{out_dir}:/work"])
        
        # Container image
        cmd.append(container_path)
        
        # ColabFold command
        cmd.append("colabfold_batch")
        cmd.extend([fasta_path, "/work"])
        
        # Add options
        cmd.extend(["--num-models", str(num_models)])
        cmd.extend(["--num-recycle", str(num_recycles)])
        
        if use_templates:
            cmd.append("--templates")
        
        if model_preset:
            # Map AlphaFold model names to ColabFold names
            colabfold_model_map = {
                "monomer": "alphafold2",
                "monomer_ptm": "alphafold2_ptm",
                "multimer": "alphafold2_multimer_v3"
            }
            mapped_preset = colabfold_model_map.get(model_preset, model_preset)
            cmd.extend(["--model-type", mapped_preset])
        
        if num_relax > 0:
            cmd.extend(["--amber", "--num-relax", str(num_relax)])
            cmd.extend(["--relax-max-iterations", str(relax_max_iteration)])
        
        if additional_options:
            cmd.extend(additional_options)
        
        return cmd

    def _build_alphafold2_native_container_command(
        self,
        fasta_path: str, # --fasta_paths
        out_dir: Path, # --output_dir
        num_relax: int = 0, # mapped to models_to_relax
        use_templates: bool = False, # if this is false, set max_template_date to 1900-01-01
        model_preset: Optional[str] = None, # --model_preset (auto default value for monomer and multimer)
        max_template_date: Optional[str] = None, # --max_template_date
        use_precomputed_msas: bool = False, # --use_precomputed_msas
        random_seed: Optional[int] = None, # --random_seed
        db_preset: str = "full_dbs", # --db_preset
        core_type: str = "gpu",
        additional_options: Optional[List[str]] = None,
        # Not exposed in af2_predict API
        num_multimer_predictions_per_model: int = 5, # --num_multimer_predictions_per_model
        benchmark: bool = False, # --benchmark
    ) -> List[str]:
        """Build command for native AlphaFold2 execution.
        This corresponds to the installation described in the AF2 repo
        (https://github.com/google-deepmind/alphafold)
        """
        # this is in principle the same as directly using run_alphafold2. The only difference is the executable name.
        # this function serve as a place holder for future difference
        cmd = self._build_alphafold2_native_python_command(
            fasta_path = fasta_path,
            out_dir = out_dir,
            num_relax = num_relax,
            use_templates = use_templates,
            model_preset = model_preset,
            max_template_date = max_template_date,
            use_precomputed_msas = use_precomputed_msas,
            random_seed = random_seed,
            db_preset = db_preset,
            core_type = core_type,
            additional_options = additional_options,
            num_multimer_predictions_per_model = num_multimer_predictions_per_model,
            benchmark = benchmark,
        )
        
        return cmd

    def _build_alphafold2_native_python_command(
        self,
        fasta_path: str, # --fasta_paths
        out_dir: Path, # --output_dir
        num_relax: int = 0, # mapped to models_to_relax
        use_templates: bool = False, # if this is false, set max_template_date to 1900-01-01
        model_preset: Optional[str] = None, # --model_preset (auto default value for monomer and multimer)
        max_template_date: Optional[str] = None, # --max_template_date
        use_precomputed_msas: bool = False, # --use_precomputed_msas
        random_seed: Optional[int] = None, # --random_seed
        db_preset: str = "full_dbs", # --db_preset
        core_type: str = "gpu",
        additional_options: Optional[List[str]] = None,
        # Not exposed in af2_predict API
        num_multimer_predictions_per_model: int = 5, # --num_multimer_predictions_per_model
        benchmark: bool = False, # --benchmark
    ) -> List[str]:
        """Build command for AlphaFold2 native Python execution.
        
        Args:
            fasta_path: Comma-separated FASTA paths passed to --fasta_paths.
            out_dir: Output directory for AlphaFold2 outputs.
            num_relax: Number of models to relax (0 none, 1 best, >1 all).
            use_templates: Whether to enable template usage.
            model_preset: AF2 model preset ('monomer_ptm' or 'multimer').
            max_template_date: Max template date when templates are used.
            use_precomputed_msas: Use precomputed MSAs (native AF2).
            random_seed: Random seed for reproducibility.
            db_preset: Database preset ('full_dbs' or 'reduced_dbs').
            core_type: Computing core type ('gpu' or 'cpu').
            additional_options: Additional CLI options to append.
            num_multimer_predictions_per_model: Multimer predictions per model.
            benchmark: Enable benchmarking mode.
        """
        afconfig = self.config_
        
        cmd = ["python", afconfig.EXECUTABLE_PATH]
        cmd.extend(["--fasta_paths", fasta_path])
        cmd.extend(["--output_dir", str(out_dir)])
        cmd.extend(["--data_dir", afconfig.DATA_DIR])
        cmd.extend(["--db_preset", db_preset])
        
        # Handle templates
        if use_templates:
            if max_template_date:
                cmd.extend(["--max_template_date", max_template_date])
            else:
                cmd.extend(["--max_template_date", "9999-12-31"])
        else:
            cmd.extend(["--max_template_date", "1900-01-01"])
        
        # Auto-detect model preset based on FASTA content if not provided
        if model_preset is None:
            # Parse FASTA to check if it's multimer
            task0_path = fasta_path.split(",")[0] # only check the first file if multiple
            fasta_sequences = parse_fasta_file(task0_path)
            if len(fasta_sequences) > 1:
                model_preset = "multimer"
            else:
                model_preset = "monomer_ptm"
        cmd.extend(["--model_preset", model_preset])
                
        # Handle relaxation based on num_relax
        if num_relax <= 0:
            cmd.extend(["--models_to_relax", "none"])
        elif num_relax == 1:
            cmd.extend(["--models_to_relax", "best"])
        else:
            _LOGGER.warning("AlphaFold2 native Python interface only supports relaxing 0, 1, or all models. Setting to relax all models.")
            cmd.extend(["--models_to_relax", "all"])
        
        # Add multimer-specific options
        if model_preset == "multimer":
            cmd.extend(["--num_multimer_predictions_per_model", str(num_multimer_predictions_per_model)])
        
        # Add other options
        if use_precomputed_msas:
            cmd.append("--use_precomputed_msas") # NOTE this currently has no effect

        if random_seed is not None:
            cmd.extend(["--random_seed", str(random_seed)])
        
        if benchmark:
            cmd.append("--benchmark")
                
        # Add all database paths required by AlphaFold using getter functions
        uniref90_path = afconfig.get_uniref90_database_path()
        if uniref90_path:
            cmd.extend(["--uniref90_database_path", uniref90_path])

        if model_preset != "multimer":
            pdb70_path = afconfig.get_pdb70_database_path()
            if pdb70_path:
                cmd.extend(["--pdb70_database_path", pdb70_path])

        mgnify_path = afconfig.get_mgnify_database_path()
        if mgnify_path:
            cmd.extend(["--mgnify_database_path", mgnify_path])
            
        uniref30_path = afconfig.get_uniref30_database_path()
        if uniref30_path:
            cmd.extend(["--uniref30_database_path", uniref30_path])
            
        bfd_path = afconfig.get_bfd_database_path()
        if bfd_path:
            cmd.extend(["--bfd_database_path", bfd_path])
            
        # small_bfd_database_path is only for reduced_dbs preset
        if db_preset == "reduced_dbs":
            small_bfd_path = afconfig.get_small_bfd_database_path()
            if small_bfd_path:
                cmd.extend(["--small_bfd_database_path", small_bfd_path])
            
        template_mmcif_dir = afconfig.get_template_mmcif_dir()
        if template_mmcif_dir:
            cmd.extend(["--template_mmcif_dir", template_mmcif_dir])
            
        # PDB seqres and Uniprot are only for multimer model preset
        if model_preset == "multimer":
            pdb_seqres_path = afconfig.get_pdb_seqres_database_path()
            if pdb_seqres_path:
                cmd.extend(["--pdb_seqres_database_path", pdb_seqres_path])
                
            uniprot_path = afconfig.get_uniprot_database_path()
            if uniprot_path:
                cmd.extend(["--uniprot_database_path", uniprot_path])

        obsolete_pdbs_path = afconfig.get_obsolete_pdbs_path()
        if obsolete_pdbs_path:
            cmd.extend(["--obsolete_pdbs_path", obsolete_pdbs_path])
        
        # Add binary paths using getter functions
        hhblits_path = afconfig.get_hhblits_binary_path()
        if hhblits_path:
            cmd.extend(["--hhblits_binary_path", hhblits_path])
            
        hhsearch_path = afconfig.get_hhsearch_binary_path()
        if hhsearch_path:
            cmd.extend(["--hhsearch_binary_path", hhsearch_path])
            
        hmmbuild_path = afconfig.get_hmmbuild_binary_path()
        if hmmbuild_path:
            cmd.extend(["--hmmbuild_binary_path", hmmbuild_path])
            
        hmmsearch_path = afconfig.get_hmmsearch_binary_path()
        if hmmsearch_path:
            cmd.extend(["--hmmsearch_binary_path", hmmsearch_path])
            
        jackhmmer_path = afconfig.get_jackhmmer_binary_path()
        if jackhmmer_path:
            cmd.extend(["--jackhmmer_binary_path", jackhmmer_path])
            
        kalign_path = afconfig.get_kalign_binary_path()
        if kalign_path:
            cmd.extend(["--kalign_binary_path", kalign_path])
        
        # Add GPU relax option
        if afconfig.USE_GPU_RELAX and core_type == "gpu":
            cmd.append("--use_gpu_relax")
        
        if additional_options:
            cmd.extend(additional_options)
        
        return cmd

    def _colabfold_format_fasta(self, sequences: Dict[str, Union[List[str], List[List[str]]]], out_dir: Union[str, Path]) -> str:
        """Convert sequences into ColabFold format FASTA file.
        
        Args:
            sequences: Mapping {seq_id: sequence} or {seq_id: [chain_seqs, ...]}.
                Accepts monomers (List[str]) or multimers (List[List[str]] per ID).
            out_dir: Output directory for the FASTA file.
            
        Returns:
            Path to the generated FASTA file (string)
        """
        fasta_path = fs.get_valid_temp_name(f"{out_dir}/input_sequences.fasta")

        if isinstance(list(sequences.values())[0], list):
            fasta_sequences = []
            sequence_ids = []

            for seq_id, multimer_seqs in sequences.items():
                combined_seq = ":".join(multimer_seqs)
                fasta_sequences.append(combined_seq)
                sequence_ids.append(seq_id)
        else:
            sequence_ids, fasta_sequences = zip(*sequences.items())

        create_fasta_from_sequences(
            fasta_sequences, 
            sequence_ids,
            output_path=fasta_path
        )

        return fasta_path

    def _native_af2_format_fasta(self, sequences: Dict[str, Union[List[str], List[List[str]]]], out_dir: Union[str, Path]) -> List[str]:
        """Create separate FASTA files for native AlphaFold input.
        
        Args:
            sequences: Dictionary mapping seq_id to per-target sequences:
                - List[str] for monomers
                - List[str] (multiple chains) for multimers
            out_dir: Output directory for FASTA files.
            
        Returns:
            List of FASTA file paths created
        """
        fasta_paths = []

        for seq_id, target_seqs in sequences.items():
            fasta_path = fs.get_valid_temp_name(f"{out_dir}/{seq_id}.fasta")
            if isinstance(target_seqs, list):
                # Create individual chain IDs for each chain in the multimer
                seq_id = [f"{seq_id}_chain_{j}" for j in range(len(target_seqs))]

            # Create FASTA file for this multimer (all chains in one file)
            create_fasta_from_sequences(
                target_seqs,
                seq_id, 
                output_path=fasta_path
            )
            fasta_paths.append(fasta_path)
        
        return fasta_paths
    
    def _validate_sequences(self, sequences: Union[List[str], List[List[str]]]) -> bool:
        """Validate that sequences are consistently either all multimers or all monomers.
        
        Args:
            sequences: Input sequences to validate
            
        Returns:
            True if sequences are multimers, False if monomers
            
        Raises:
            ValueError: If sequences are mixed or invalid format
        """
        if not sequences:
            _LOGGER.error("Empty sequences list provided")
            raise ValueError("Empty sequences list provided")

        if not isinstance(sequences, list):
            _LOGGER.error("Invalid sequences format: Expected list")
            raise ValueError("Invalid sequences format: Expected list")

        # Check first sequence to determine expected type
        first_is_list = isinstance(sequences[0], list)
        
        # Validate all sequences have consistent type
        for i, seq in enumerate(sequences):
            is_list = isinstance(seq, list)
            if is_list != first_is_list:
                _LOGGER.error(f"Inconsistent sequence format: sequence {i} is {'list' if is_list else 'string'} "
                             f"but expected {'list' if first_is_list else 'string'} based on first sequence")
                raise ValueError(f"Mixed sequence formats not supported: sequence {i} format mismatch")
                
            if is_list and not seq:
                _LOGGER.error(f"Empty multimer sequence at index {i}")
                raise ValueError(f"Empty multimer sequence at index {i}")
                
        return first_is_list

    def _parse_comprehensive_results(self, work_dir: Union[str, Path], seq_id: str, parser: PDBParser) -> Dict:
        """Parse comprehensive AlphaFold results for a sequence including all models and scores.
        
        Supports both ColabFold and native AF2 output formats.
        
        Args:
            work_dir: Directory containing AlphaFold output files
            seq_id: Sequence identifier (e.g., 'seq_0')
            parser: PDB parser instance
            
        Returns:
            Dict containing comprehensive results with all models and scores
        """
        work_dir = Path(work_dir)
        
        # Check if this is native AF2 format (has ranking_debug.json and .pkl files)
        seq_subdir = work_dir / seq_id
        ranking_debug_file = seq_subdir / "ranking_debug.json"
        is_native_af2 = ranking_debug_file.exists() and list(seq_subdir.glob("result_model_*.pkl"))
        
        if is_native_af2:
            return self._parse_native_af2_results(seq_subdir, seq_id, parser)
        else:
            return self._parse_colabfold_results(work_dir, seq_id, parser)
    
    def _parse_native_af2_results(self, seq_dir: Path, seq_id: str, parser: PDBParser) -> Dict:
        """Parse native AlphaFold2 results format.
        
        Args:
            seq_dir: Directory containing sequence results (work_dir/seq_id/)
            seq_id: Sequence identifier
            parser: PDB parser instance
            
        Returns:
            Dict containing comprehensive results
        """
        results = {}
        
        # Load ranking information
        ranking_file = seq_dir / "ranking_debug.json"
        with open(ranking_file, 'r') as f:
            ranking_data = json.load(f)
        model_order = ranking_data['order']  # Best to worst
        
        # Find all available models
        pkl_files = list(seq_dir.glob("result_model_*.pkl"))
        ranked_pdb_files = list(seq_dir.glob("ranked_*.pdb"))
        
        if not pkl_files or not ranked_pdb_files:
            _LOGGER.error(f"No pkl files or ranked PDB files found for sequence {seq_id}")
            raise ValueError(f"No pkl files or ranked PDB files found for sequence {seq_id}")
        
        # Parse structures and scores for each model
        model_data = {}
        
        # Load pLDDT scores from pkl files
        for pkl_file in pkl_files:
            # Extract model name (e.g., model_1_ptm_pred_0)
            match = re.search(r'result_(model_\d+.+)\.pkl', pkl_file.name)
            if match:
                model_name = match.group(1)
                pkl_data = load_obj(str(pkl_file))
                plddt_scores = pkl_data.get('plddt', [])
                model_data[model_name] = {
                    'plddt': plddt_scores.tolist() if hasattr(plddt_scores, 'tolist') else list(plddt_scores)
                }
        
        # Load structures from ranked PDB files
        for i, ranked_pdb in enumerate(sorted(ranked_pdb_files)):
            try:
                structure = parser.get_structure(str(ranked_pdb))
                rank = i  # 0-based ranking
                
                if rank < len(model_order):
                    model_name = model_order[rank]
                    if model_name in model_data:
                        model_data[model_name]['structure'] = structure
                        model_data[model_name]['rank'] = rank
            except Exception as e:
                _LOGGER.warning(f"Could not parse structure file {ranked_pdb}: {e}")
        
        if not model_data:
            _LOGGER.error(f"No valid model data found for sequence {seq_id}")
            raise ValueError(f"No valid model data found for sequence {seq_id}")
        
        # Build comprehensive results
        # Best model is first in order
        best_model_name = model_order[0]
        if best_model_name in model_data and 'structure' in model_data[best_model_name]:
            # Extract model number (e.g., 5 from model_5_ptm_pred_0)
            model_num_match = re.search(r'model_(\d+)', best_model_name)
            best_model_num = int(model_num_match.group(1)) if model_num_match else 1
            
            results['best_model'] = model_data[best_model_name]['structure']
            results['best_model_plddt'] = model_data[best_model_name]['plddt']
            results['best_model_index'] = best_model_num
        
        # Add all individual models
        for model_name, data in model_data.items():
            if 'structure' in data:
                model_num_match = re.search(r'model_(\d+)', model_name)
                if model_num_match:
                    model_num = int(model_num_match.group(1))
                    results[f'model_{model_num}'] = data['structure']
                    results[f'model_{model_num}_plddt'] = data['plddt']
        
        return results
    
    def _parse_colabfold_results(self, work_dir: Path, seq_id: str, parser: PDBParser) -> Dict:
        """Parse ColabFold results format.
        
        Args:
            work_dir: Directory containing AlphaFold output files
            seq_id: Sequence identifier (e.g., 'seq_0')
            parser: PDB parser instance
            
        Returns:
            Dict containing comprehensive results
        """
        # Find all PDB files for this sequence
        pdb_files = list(work_dir.glob(f"{seq_id}_*.pdb"))
        
        # Find all JSON score files for this sequence  
        json_files = list(work_dir.glob(f"{seq_id}_scores_*.json"))
        
        if not pdb_files:
            _LOGGER.warning(f"No PDB files found for sequence {seq_id}")
            return {}
        
        # Parse scores from JSON files
        model_scores = {}
        for json_file in json_files:
            # Extract model number and rank from filename
            # Format: seq_0_scores_rank_001_alphafold2_ptm_model_5_seed_000.json
            match = re.search(r'rank_(\d+).*model_(\d+)', json_file.name)
            if match:
                rank = int(match.group(1))
                model_num = int(match.group(2))
                
                with open(json_file, 'r') as f:
                    score_data = json.load(f)
                model_scores[rank] = {
                    'model_num': model_num,
                    'plddt': score_data.get('plddt', [])
                }
        
        # Parse structures from PDB files
        model_structures = {}
        for pdb_file in pdb_files:
            # Extract rank and model number from filename
            match = re.search(r'rank_(\d+).*model_(\d+)', pdb_file.name)
            if match:
                rank = int(match.group(1))
                model_num = int(match.group(2))
                
                structure = parser.get_structure(str(pdb_file))
                model_structures[rank] = {
                    'model_num': model_num,
                    'structure': structure
                }
        
        if not model_structures:
            _LOGGER.error(f"No valid structures found for sequence {seq_id}")
            return {}
        
        # Build comprehensive results
        results = {}
        
        # Get best model (rank 1)
        if 1 in model_structures and 1 in model_scores:
            best_model_num = model_structures[1]['model_num']
            results['best_model'] = model_structures[1]['structure']
            results['best_model_plddt'] = model_scores[1]['plddt']
            results['best_model_index'] = best_model_num
        
        # Add all individual models
        for rank in sorted(model_structures.keys()):
            if rank in model_scores:
                model_num = model_structures[rank]['model_num']
                results[f'model_{model_num}'] = model_structures[rank]['structure']
                results[f'model_{model_num}_plddt'] = model_scores[rank]['plddt']
        
        return results

    def _create_sequence_mapper(self, sequences: Union[List[str], List[List[str]]], work_dir: Path) -> Dict[str, Union[str, List[str]]]:
        """Creates a mapper from unique sequence IDs to sequences.

        Ensures that the generated sequence IDs do not conflict with existing
        FASTA files in the working directory.

        Args:
            sequences: A list of sequences or a list of list of sequences (for multimers).
            work_dir: The directory where output files will be stored.

        Returns:
            A dictionary mapping unique sequence IDs to their corresponding sequences.
        """
        work_dir = Path(work_dir)
        sequences_mapper = {}
        i = 0
        for seq in sequences:
            while True:
                seq_id = f"seq_{i}"
                # Check for existing FASTA file for this seq_id
                potential_fasta_path = work_dir / f"{seq_id}.fasta"
                if not potential_fasta_path.exists():
                    sequences_mapper[seq_id] = seq
                    i += 1
                    break
                i += 1
                if i > 99999: # safety
                    _LOGGER.error("Loop exceeded maximum iterations. Failed to find unique sequence ID.")
                    raise RuntimeError("Failed to generate unique sequence IDs after 9999 attempts.")
        return sequences_mapper
