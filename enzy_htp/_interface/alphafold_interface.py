"""Defines an AlphafoldInterface class that serves as a bridge for enzy_htp to utilize AlphaFold2 software.

Author: Gemini, Claude Code
Date: 2025-08-16
"""
from __future__ import annotations
import tempfile
import os
import re
from typing import Union, List, Optional, Dict, Tuple
from pathlib import Path
from dataclasses import dataclass

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
        sequences: List[str], 
        out_dir: Union[str, Path, None] = None,
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
        **kwargs
    ) -> Dict[str, Structure]:
        """Science API for AlphaFold2 structure prediction.

        Args:
            sequences: List of amino acid sequences to predict
            out_dir: Output directory for results. If None, creates temporary directory.
            cluster_job_config: Configuration for cluster job submission
            array_size: Number of jobs to run simultaneously (for cluster submission)
            job_check_period: Time cycle for job state checking (seconds)
            seq_per_job: Number of sequences per job (for array execution)
            num_models: Number of models to generate
            num_recycles: Number of recycling iterations
            num_relax: Number of structures to relax
            relax_max_iteration: Maximum relaxation iterations
            use_templates: Whether to use templates
            max_template_date: Maximum template date
            model_preset: Model preset configuration
            db_preset: Database preset
            additional_options: Additional command-line options
            **kwargs: Additional arguments
            
        Returns:
            Dict mapping sequence identifiers to Structure objects
        """
        temp_paths = []
        stru_parser = PDBParser()
        
        try:
            # Create FASTA file using helper function
            sequence_ids = [f"seq_{i}" for i in range(len(sequences))]
            fasta_path = create_fasta_from_sequences(
                sequences, 
                sequence_ids, 
                delete_on_close=False
            )
            temp_paths.append(fasta_path)

            # Create output directory
            if out_dir is None:
                out_dir = Path(tempfile.mkdtemp(prefix='alphafold2_'))
                temp_paths.append(str(out_dir))
            else:
                out_dir = Path(out_dir)
            out_dir.mkdir(exist_ok=True)

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
                    additional_options=additional_options
                )

            # Parse results into Structure objects
            # Map sequence IDs back to full sequences for return values
            fasta_sequences = parse_fasta_file(fasta_path)
            seq_id_to_sequence = {seq_id: seq for seq_id, seq in fasta_sequences}
            
            structures = {}
            for seq_id, pdb_path in result_files.items():
                structure = stru_parser.get_structure(pdb_path)
                # Use full sequence as key instead of sequence ID
                full_sequence = seq_id_to_sequence[seq_id]
                structures[full_sequence] = structure

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
                model_preset, additional_options
            )
        elif config.INSTALL_TYPE == "alphafold2_native_container":
            cmd = self._build_alphafold2_native_container_command(
                fasta_path, out_dir, max_template_date, 
                model_preset, db_preset, additional_options
            )
        elif config.INSTALL_TYPE == "alphafold2_native_python":
            cmd = self._build_alphafold2_native_python_command(
                fasta_path, out_dir, num_models, num_recycles,
                num_relax, relax_max_iteration, use_templates,
                model_preset, additional_options
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
            
            # Create FASTA for this job
            job_fasta_path = create_fasta_from_sequences(
                [seq for _, seq in job_sequences],
                job_seq_ids,
                output_path=out_dir / f"job_{i//seq_per_job}.fasta"
            )
            
            # Build command for this job
            config = self.config_
            
            if config.INSTALL_TYPE == "colabfold_container":
                cmd = self._build_colabfold_container_command(
                    job_fasta_path, out_dir, num_models, num_recycles, 
                    num_relax, relax_max_iteration, use_templates, 
                    model_preset, additional_options
                )
            elif config.INSTALL_TYPE == "alphafold2_native_container":
                cmd = self._build_alphafold2_native_container_command(
                    job_fasta_path, out_dir, max_template_date, 
                    model_preset, db_preset, additional_options
                )
            elif config.INSTALL_TYPE == "alphafold2_native_python":
                cmd = self._build_alphafold2_native_python_command(
                    job_fasta_path, out_dir, num_models, num_recycles,
                    num_relax, relax_max_iteration, use_templates,
                    model_preset, additional_options
                )
            else:
                raise ValueError(f"Unsupported install type: {config.INSTALL_TYPE}")
            
            # Create ClusterJob
            if isinstance(cluster_job_config, dict):
                job_config = ClusterJobConfig.from_dict(cluster_job_config)
            else:
                job_config = cluster_job_config
            
            # Use ClusterJob.config_job to create the job properly
            if not job_config.has_cluster():
                raise ValueError("cluster_job_config must specify a cluster for job execution")
            
            if not job_config.has_res_keywords():
                _LOGGER.warning("No resource keywords specified in cluster job config, using empty dict")
                job_config.res_keywords = {}
            
            job = ClusterJob.config_job(
                commands=' '.join(cmd),
                cluster=job_config.cluster,
                env_settings=[],
                res_keywords=job_config.res_keywords,
                sub_dir=str(out_dir)
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
        additional_options: Optional[List[str]]
    ) -> List[str]:
        """Build command for ColabFold container execution."""
        config = self.config_
        
        # Expand user paths
        container_path = os.path.expanduser(config.CONTAINER_PATH)
        
        cmd = [config.CONTAINER_TYPE, "run"]
        
        # Add GPU support if available
        if config.CONTAINER_TYPE in ["docker", "apptainer", "singularity"]:
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
        additional_options: Optional[List[str]]
    ) -> List[str]:
        """Build command for native AlphaFold2 execution."""
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
        
        if hasattr(config, 'USE_GPU_RELAX') and config.USE_GPU_RELAX:
            cmd.append("--use_gpu_relax")
        
        if additional_options:
            cmd.extend(additional_options)
        
        return cmd

    def _build_alphafold2_native_python_command(
        self,
        fasta_path: str,
        out_dir: Path,
        num_models: int,
        num_recycles: int,
        num_relax: int,
        relax_max_iteration: int,
        use_templates: bool,
        model_preset: str,
        additional_options: Optional[List[str]]
    ) -> List[str]:
        """Build command for AlphaFold2 Python execution."""
        config = self.config_
        
        cmd = [config.EXECUTABLE_PATH]
        cmd.extend([fasta_path, str(out_dir)])
        
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

    def _get_rank(self, filename):
        match = re.search(r'rank_(\d+)', filename)
        return int(match.group(1)) if match else 999  # Lower rank number is better
