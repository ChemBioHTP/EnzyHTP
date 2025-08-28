"""Defines an AlphafoldInterface class that serves as a bridge for enzy_htp to utilize AlphaFold2 software.

Author: Gemini, Claude Code
Date: 2025-08-16
"""
from __future__ import annotations
import tempfile
import os
import re
import copy
from typing import Union, List, Optional, Dict, Tuple
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
    
    def get_expected_output_files(self) -> Dict[str, str]:
        """Get expected output PDB file paths for each sequence.
        
        Returns:
            Dict mapping sequence IDs to expected PDB file paths
        """
        # Get all available files first
        filename_to_path = self.interface._find_output_files_map(self.output_dir)
        
        # Filter and select best files for this job's sequences
        return self.interface._select_best_files_for_sequences(filename_to_path, self.sequence_ids)


class AlphafoldInterface(BaseInterface):
    """Class that provides a direct interface for enzy_htp to utilize AlphaFold2 software."""

    def __init__(self, parent, config: AlphafoldConfig = None) -> None:
        """Simplistic constructor that optionally takes an AlphafoldConfig object as its only argument.
        Calls parent class."""
        super().__init__(parent, config, AlphafoldConfig)

    def af2_predict(
        self, 
        sequences: Union[List[str], List[List[str]]], 
        out_dir: Union[str, Path, None] = None,
        # local run related
        non_armer_core_type: str = "gpu",
        # cluster job related
        cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None,
        array_size: int = 0,
        job_check_period: int = 30,
        seq_per_job: int = 1,
        # AlphaFold2 specific
        num_models: int = 5,
        num_recycles: int = 3,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = None,
        db_preset: str = "reduced_dbs",
        additional_options: Optional[List[str]] = None,
        use_precomputed_msas: bool = False,
        random_seed: Optional[int] = None,
        **kwargs
    ) -> Dict[str, Structure]:
        """Science API for AlphaFold2 structure prediction.

        Args:
            sequences: Sequences to predict. Can be:
                      - List[str]: Single sequences [seq1, seq2, seq3, ...]
                      - List[List[str]]: Multimers [[seq1_chain_A, seq1_chain_B], [seq2_chain_A, seq2_chain_B], ...]
            out_dir: Output directory for results. If None, creates temporary directory.
            cluster_job_config: Configuration for cluster job submission
            core_type: Type of computing core ('gpu' or 'cpu') for cluster jobs
            array_size: Number of jobs to run simultaneously (for cluster submission)
            job_check_period: Time cycle for job state checking (seconds)
            seq_per_job: Number of sequences per job (for array execution)
            num_models: Number of models to generate
            num_recycles: Number of recycling iterations
            num_relax: Number of top ranked structures to relax
            relax_max_iteration: Maximum relaxation iterations
            use_templates: Whether to use templates
            max_template_date: Maximum template date
            model_preset: Model preset configuration
            db_preset: Database preset
            additional_options: Additional command-line options
            use_precomputed_msas: Whether to use precomputed MSAs 
                (NOTE: when native af2 distribution is used. HHSearch will still be run. Details: https://github.com/google-deepmind/alphafold/issues/469)
            random_seed: Random seed for reproducibility
            
        Returns:
            Dict mapping sequence identifiers to Structure objects
        """
        temp_paths = []
        stru_parser = PDBParser()
        
        try:
            # Create output directory first
            if out_dir is None:
                scratch_dir = eh_config.system.SCRATCH_DIR
                out_dir = fs.get_valid_temp_name(f"{scratch_dir}/alphafold2")
                temp_paths.append(out_dir)
            fs.safe_mkdir(out_dir)

            # Handle different sequence input formats
            fasta_sequences, sequence_ids = self._prepare_sequences_and_ids(sequences)
            fasta_path = fs.get_valid_temp_name(f"{out_dir}/input_sequences.fasta")
            create_fasta_from_sequences(
                fasta_sequences, 
                sequence_ids,
                output_path=fasta_path
            )
            temp_paths.append(fasta_path)

            if cluster_job_config:
                # Run on cluster using array jobs
                result_eggs = self.make_job(
                    fasta_path=fasta_path,
                    out_dir=out_dir,
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
                    seq_per_job=seq_per_job
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
                
                # Collect results from all eggs
                result_files = {}
                for egg in result_eggs:
                    egg_files = egg.get_expected_output_files()
                    result_files.update(egg_files)
            else:
                # Run locally
                result_files = self.run(
                    fasta_path=fasta_path,
                    out_dir=out_dir,
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
                    non_armer_core_type=non_armer_core_type
                )

            # Parse results into Structure objects
            # Map sequence IDs back to original input format
            structures = {}
            for seq_id, pdb_path in result_files.items():
                structure = stru_parser.get_structure(pdb_path)
                # Map back to original sequence format
                original_seq = self._map_seq_id_to_original(seq_id, sequences)
                structures[original_seq] = structure

            return structures

        finally:
            # Clean up temporary files
            fs.clean_temp_file_n_dir(temp_paths)

    def run(
        self,
        fasta_path: str,
        out_dir: Union[str, Path],
        num_models: int = 5,
        num_recycles: int = 3,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = "alphafold2_ptm",
        db_preset: str = "reduced_dbs",
        additional_options: Optional[List[str]] = None,
        use_precomputed_msas: bool = False,
        random_seed: Optional[int] = None,
        non_armer_core_type: str = "gpu",
    ) -> Dict[str, str]:
        """Execute AlphaFold2 prediction locally.
        
        Args:
            fasta_path: Path to input FASTA file
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
            Dict mapping sequence IDs to output PDB file paths
        """
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)
        
        config = self.config_
        
        # build command
        if config.INSTALL_TYPE == "colabfold_container":
            cmd = self._build_colabfold_container_command(
                fasta_path, out_dir, num_models, num_recycles, 
                num_relax, relax_max_iteration, use_templates, 
                model_preset, additional_options, non_armer_core_type
            )
        elif config.INSTALL_TYPE == "alphafold2_native_container":
            cmd = self._build_alphafold2_native_container_command(
                fasta_path, out_dir, max_template_date, 
                model_preset, db_preset, additional_options, non_armer_core_type
            )
        elif config.INSTALL_TYPE == "alphafold2_native_python":
            cmd = self._build_alphafold2_native_python_command(
                fasta_path, out_dir, num_relax, use_templates,
                model_preset, additional_options, 
                max_template_date=max_template_date,
                db_preset=db_preset,
                use_precomputed_msas=use_precomputed_msas,
                random_seed=random_seed,
                core_type=non_armer_core_type
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
        
        # Get all files and select best ones for sequences
        filename_to_path = self._find_output_files_map(out_dir)
        fasta_sequences = parse_fasta_file(fasta_path)
        sequence_ids = [seq_id for seq_id, _ in fasta_sequences]
        return self._select_best_files_for_sequences(filename_to_path, sequence_ids)

    def make_job(
        self,
        fasta_path: str,
        out_dir: Union[str, Path],
        num_models: int = 5,
        num_recycles: int = 3,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = "alphafold2_ptm",
        db_preset: str = "reduced_dbs",
        additional_options: Optional[List[str]] = None,
        use_precomputed_msas: bool = False,
        random_seed: Optional[int] = None,
        cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None,
        seq_per_job: int = 1,
    ) -> List[AlphaFold2ResultEgg]:
        """Create cluster jobs for running AlphaFold2 with array support.
        
        Args:
            fasta_path: Path to input FASTA file
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
            cluster_job_config: Configuration for cluster job
            seq_per_job: Number of sequences per job
            
        Returns:
            List of AlphaFold2ResultEgg objects, one per job
        """
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)
        
        # Parse FASTA to get sequences
        fasta_sequences = parse_fasta_file(fasta_path)
        
        # Create individual jobs for array execution
        result_eggs = []
        
        for i in range(0, len(fasta_sequences), seq_per_job):
            job_sequences = fasta_sequences[i:i+seq_per_job]
            job_seq_ids = [seq_id for seq_id, _ in job_sequences]
            job_seq = [seq for _, seq in job_sequences]
            
            # Create FASTA for this job
            job_fasta_path = create_fasta_from_sequences(
                job_seq,
                job_seq_ids,
                output_path=out_dir / f"job_{i//seq_per_job}.fasta"
            )
            
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
            config = self.config_
            
            if config.INSTALL_TYPE == "colabfold_container":
                cmd = self._build_colabfold_container_command(
                    job_fasta_path, out_dir, num_models, num_recycles, 
                    num_relax, relax_max_iteration, use_templates, 
                    model_preset, additional_options, core_type
                )
            elif config.INSTALL_TYPE == "alphafold2_native_container":
                cmd = self._build_alphafold2_native_container_command(
                    job_fasta_path, out_dir, max_template_date, 
                    model_preset, db_preset, additional_options, core_type
                )
            elif config.INSTALL_TYPE == "alphafold2_native_python":
                cmd = self._build_alphafold2_native_python_command(
                    job_fasta_path, out_dir, num_relax, use_templates,
                    model_preset, additional_options,
                    max_template_date=max_template_date,
                    db_preset=db_preset,
                    use_precomputed_msas=use_precomputed_msas,
                    random_seed=random_seed,
                    core_type=core_type
                )
            else:
                raise ValueError(f"Unsupported install type: {config.INSTALL_TYPE}")
            
            # Handle default res_keywords similar to amber_interface
            res_keywords_update = job_config.res_keywords if job_config.has_res_keywords() else {}
            default_res_keywords = self.config_.get_default_af2_cluster_job_res_keywords(core_type)
            job_config.res_keywords = default_res_keywords | res_keywords_update
            
            # Set up environment settings based on cluster and core type
            cluster = job_config.cluster
            env_settings = cluster.AF2_ENV[core_type.upper()]
            
            # Create submission script path
            sub_script_path = out_dir / f"submit_alphafold_{i//seq_per_job}.cmd"
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
            non_armer_core_type: Computing core type ('gpu' or 'cpu')
        """
        config = self.config_
        
        # Expand user paths
        container_path = os.path.expanduser(config.CONTAINER_PATH)
        
        cmd = [config.CONTAINER_TYPE, "run"]
        
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
        fasta_path: str,
        out_dir: Path,
        max_template_date: Optional[str],
        model_preset: str,
        db_preset: str,
        additional_options: Optional[List[str]],
        core_type: str = "gpu"
    ) -> List[str]:
        """Build command for native AlphaFold2 execution.
        
        Args:
            non_armer_core_type: Computing core type ('gpu' or 'cpu')
        """
        config = self.config_
        
        cmd = ["python", config.EXECUTABLE_PATH]
        cmd.extend(["--fasta_paths", fasta_path])
        cmd.extend(["--output_dir", str(out_dir)])
        cmd.extend(["--data_dir", config.DATA_DIR])
        cmd.extend(["--model_preset", model_preset])
        cmd.extend(["--db_preset", db_preset])
        
        if max_template_date:
            cmd.extend(["--max_template_date", max_template_date])
        
        # Add database paths (these would typically be in config)
        if hasattr(config, 'UNIREF90_DATABASE_PATH') and config.UNIREF90_DATABASE_PATH:
            cmd.extend(["--uniref90_database_path", config.UNIREF90_DATABASE_PATH])
        
        if hasattr(config, 'MGNIFY_DATABASE_PATH') and config.MGNIFY_DATABASE_PATH:
            cmd.extend(["--mgnify_database_path", config.MGNIFY_DATABASE_PATH])
        
        if hasattr(config, 'TEMPLATE_MMCIF_DIR') and config.TEMPLATE_MMCIF_DIR:
            cmd.extend(["--template_mmcif_dir", config.TEMPLATE_MMCIF_DIR])
        
        if hasattr(config, 'USE_GPU_RELAX') and config.USE_GPU_RELAX and core_type == "gpu":
            cmd.append("--use_gpu_relax")
        
        if additional_options:
            cmd.extend(additional_options)
        
        return cmd

    def _build_alphafold2_native_python_command(
        self,
        fasta_path: str, # --fasta_paths
        out_dir: Path, # --output_dir
        num_relax: int = 0, # mapped to models_to_relax
        use_templates: bool = False, # if this is false, set max_template_date to 1900-01-01
        model_preset: Optional[str] = None, # --model_preset (auto default value for monomer and multimer)
        additional_options: Optional[List[str]] = None,
        max_template_date: Optional[str] = None, # --max_template_date
        num_multimer_predictions_per_model: int = 5, # --num_multimer_predictions_per_model
        use_precomputed_msas: bool = False, # --use_precomputed_msas
        random_seed: Optional[int] = None, # --random_seed
        db_preset: str = "full_dbs", # --db_preset
        benchmark: bool = False, # --benchmark
        core_type: str = "gpu",
    ) -> List[str]:
        """Build command for AlphaFold2 native Python execution.
        
        Args:
            non_armer_core_type: Computing core type ('gpu' or 'cpu')
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
            fasta_sequences = parse_fasta_file(fasta_path)
            if len(fasta_sequences) > 1 or any(':' in seq for _, seq in fasta_sequences): # BUG native AF2 should not take : containing sequence
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
                
        # Add all database paths required by AlphaFold
        if afconfig.UNIREF90_DATABASE_PATH:
            cmd.extend(["--uniref90_database_path", afconfig.UNIREF90_DATABASE_PATH])

        if afconfig.PDB70_DATABASE_PATH:
            cmd.extend(["--pdb70_database_path", afconfig.PDB70_DATABASE_PATH])

        if afconfig.MGNIFY_DATABASE_PATH:
            cmd.extend(["--mgnify_database_path", afconfig.MGNIFY_DATABASE_PATH])
            
        if afconfig.UNIREF30_DATABASE_PATH:
            cmd.extend(["--uniref30_database_path", afconfig.UNIREF30_DATABASE_PATH])
            
        if afconfig.BFD_DATABASE_PATH:
            cmd.extend(["--bfd_database_path", afconfig.BFD_DATABASE_PATH])
            
        # small_bfd_database_path is only for reduced_dbs preset
        if db_preset == "reduced_dbs" and afconfig.SMALL_BFD_DATABASE_PATH:
            cmd.extend(["--small_bfd_database_path", afconfig.SMALL_BFD_DATABASE_PATH])
            
        if afconfig.TEMPLATE_MMCIF_DIR:
            cmd.extend(["--template_mmcif_dir", afconfig.TEMPLATE_MMCIF_DIR])
            
        # PDB seqres and Uniprot are only for multimer model preset
        if model_preset == "multimer":
            if afconfig.PDB_SEQRES_DATABASE_PATH:
                cmd.extend(["--pdb_seqres_database_path", afconfig.PDB_SEQRES_DATABASE_PATH])
                
            if afconfig.UNIPROT_DATABASE_PATH:
                cmd.extend(["--uniprot_database_path", afconfig.UNIPROT_DATABASE_PATH])

        if afconfig.OBSOLETE_PDBS_PATH:
            cmd.extend(["--obsolete_pdbs_path", afconfig.OBSOLETE_PDBS_PATH])
        
        # Add binary paths
        if afconfig.HHBLITS_BINARY_PATH:
            cmd.extend(["--hhblits_binary_path", afconfig.HHBLITS_BINARY_PATH])
            
        if afconfig.HHSEARCH_BINARY_PATH:
            cmd.extend(["--hhsearch_binary_path", afconfig.HHSEARCH_BINARY_PATH])
            
        if afconfig.HMMBUILD_BINARY_PATH:
            cmd.extend(["--hmmbuild_binary_path", afconfig.HMMBUILD_BINARY_PATH])
            
        if afconfig.HMMSEARCH_BINARY_PATH:
            cmd.extend(["--hmmsearch_binary_path", afconfig.HMMSEARCH_BINARY_PATH])
            
        if afconfig.JACKHMMER_BINARY_PATH:
            cmd.extend(["--jackhmmer_binary_path", afconfig.JACKHMMER_BINARY_PATH])
            
        if afconfig.KALIGN_BINARY_PATH:
            cmd.extend(["--kalign_binary_path", afconfig.KALIGN_BINARY_PATH])
        
        # Add GPU relax option
        if afconfig.USE_GPU_RELAX and core_type == "gpu":
            cmd.append("--use_gpu_relax")
        
        if additional_options:
            cmd.extend(additional_options)
        
        return cmd

    def _prepare_sequences_and_ids(self, sequences: Union[List[str], List[List[str]]]) -> Tuple[List[str], List[str]]:
        """Prepare sequences and IDs for FASTA creation based on install type.
        
        Args:
            sequences: Either single sequences or multimer sequences
            
        Returns:
            Tuple of (fasta_sequences, sequence_ids)
        """
        config = self.config_
        
        if isinstance(sequences[0], list):
            # Multimer format: [[seq1_chain_A, seq1_chain_B], [seq2_chain_A, seq2_chain_B], ...]
            fasta_sequences = []
            sequence_ids = []
            
            for i, multimer_seqs in enumerate(sequences):
                seq_id = f"seq_{i}"
                sequence_ids.append(seq_id)
                
                if config.INSTALL_TYPE == "colabfold_container":
                    # ColabFold: use ':' to separate chains in same sequence
                    combined_seq = ":".join(multimer_seqs)
                    fasta_sequences.append(combined_seq)
                else:
                    # Native AlphaFold: put multiple sequences in same FASTA file
                    # Only add the first sequence here, as native AF handles multiple sequences in file
                    fasta_sequences.append(":".join(multimer_seqs))
        else:
            # Single sequence format: [seq1, seq2, seq3, ...]
            fasta_sequences = list(sequences)
            sequence_ids = [f"seq_{i}" for i in range(len(sequences))]
            
        return fasta_sequences, sequence_ids
    
    def _map_seq_id_to_original(self, seq_id: str, original_sequences: Union[List[str], List[List[str]]]) -> Union[str, Tuple[str, ...]]:
        """Map sequence ID back to original input format.
        
        Args:
            seq_id: Sequence ID like 'seq_0', 'seq_1', etc.
            original_sequences: Original input sequences
            
        Returns:
            Original sequence or tuple of sequences for multimers
        """
        try:
            index = int(seq_id.split('_')[1])
            original_seq = original_sequences[index]
            
            if isinstance(original_seq, list):
                return tuple(original_seq)  # Return as tuple for hashable key
            else:
                return original_seq
        except (IndexError, ValueError):
            _LOGGER.warning(f"Could not map sequence ID {seq_id} back to original sequence")
            return seq_id

    def _find_output_files_map(self, out_dir: Path) -> Dict[str, str]:
        """Find output PDB files in the output directory and return filename-to-path mapping.
        
        Args:
            out_dir: Output directory to search
            
        Returns:
            Dict mapping filenames to PDB file paths
        """
        filename_to_path = {}
        
        # Look for PDB files
        for pdb_file in out_dir.glob("**/*.pdb"):
            filename = pdb_file.stem
            filename_to_path[filename] = str(pdb_file)
        
        return filename_to_path
    
    def _select_best_files_for_sequences(self, filename_to_path: Dict[str, str], sequence_ids: List[str]) -> Dict[str, str]:
        """Filter and select the best files for the given sequences.
        
        Args:
            filename_to_path: Mapping of filename to file path
            sequence_ids: List of sequence identifiers to find files for
            
        Returns:
            Dict mapping sequence IDs to the best PDB file path for each sequence
        """
        result_files = {}
        
        for seq_id in sequence_ids:
            # Find all files for this sequence
            seq_files = {}
            
            for filename, filepath in filename_to_path.items():
                # Check if this file belongs to the current sequence
                if filename.startswith(f"{seq_id}_"):
                    seq_files[filename] = filepath
            
            if seq_files:
                # Select the best file: prefer relaxed over unrelaxed, and higher ranked models
                best_filename = self._select_best_file(seq_files.keys())
                result_files[seq_id] = seq_files[best_filename]
            else:
                _LOGGER.warning(f"No files found for sequence ID {seq_id}. "
                "AlphaFold may have error on sequence.")

        return result_files
    
    def _select_best_file(self, filenames: List[str]) -> str:
        """Select the best filename from a list based on ranking and relaxation status.
        
        Args:
            filenames: List of filenames to choose from
            
        Returns:
            The best filename
        """
        if not filenames:
            _LOGGER.error("No filenames provided to select best file")
            raise ValueError("No filenames provided")
            
        if len(filenames) == 1:
            return list(filenames)[0]
        
        relaxed_files = []
        unrelaxed_files = []
        
        for filename in filenames:
            if "unrelaxed" in filename:
                unrelaxed_files.append(filename)
            else:
                relaxed_files.append(filename)
        
        # Prefer relaxed files
        candidates = relaxed_files if relaxed_files else unrelaxed_files
        
        # If we have multiple candidates, prefer higher ranked models
        # Look for patterns like rank_001, rank_002, etc.
        return min(candidates, key=self._get_rank)

    def _get_rank(self, filename: str) -> int:
        """Extract rank number from filename for sorting purposes.
        
        Args:
            filename: Filename containing rank information (e.g., 'seq_0_rank_001_model_1.pdb')
            
        Returns:
            Rank number if found, otherwise 999 (for consistent sorting behavior)
        """
        match = re.search(r'rank_(\d+)', filename)
        return int(match.group(1)) if match else 999  # Lower rank number is better
