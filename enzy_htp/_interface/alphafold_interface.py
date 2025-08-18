"""Defines an AlphafoldInterface class that serves as a bridge for enzy_htp to utilize Alphafold software.

Author: Gemini
Date: 2025-08-16
"""
from __future__ import annotations
import tempfile
import subprocess
import os
from typing import Union, List, Optional, Dict
from pathlib import Path

from .base_interface import BaseInterface
from enzy_htp.structure import Structure
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp._config.alphafold_config import AlphafoldConfig
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig
from enzy_htp.core.logger import _LOGGER

class AlphafoldInterface(BaseInterface):
    """Class that provides a direct inteface for enzy_htp to utilize Alphafold software.
    """

    def __init__(self, parent, config: AlphafoldConfig = None) -> None:
        """Simplistic constructor that optionally takes an AlphafoldConfig object as its only argument.
        Calls parent class."""
        super().__init__(parent, config, AlphafoldConfig)

    def af2_predict(self, sequences: List[str], cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None, **kwargs) -> Dict[str, Structure]:
        """Wrapper for running AlphaFold2 prediction.

        This function is a wrapper for the `run` method. It handles the submission of the job to a cluster if
        `cluster_job_config` is provided.
        
        Args:
            sequences: List of amino acid sequences to predict
            cluster_job_config: Configuration for cluster job submission
            **kwargs: Additional arguments passed to run() or make_job()
            
        Returns:
            Dict mapping sequence identifiers to Structure objects
        """
        parser = PDBParser()
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as fasta_file:
            for i, seq in enumerate(sequences):
                fasta_file.write(f">seq_{i}\n{seq}\n")
            fasta_path = fasta_file.name

        try:
            # Create output directory
            out_dir = Path(kwargs.get('out_dir', tempfile.mkdtemp(prefix='alphafold_')))
            out_dir.mkdir(exist_ok=True)

            if cluster_job_config:
                # Run on cluster
                job = self.make_job(fasta_path, out_dir, **kwargs)
                job.submit()
                job.wait_to_finish()
                if job.get_status() != "completed":
                    raise RuntimeError(f"AlphaFold job failed with status: {job.get_status()}")
                result_files = self._find_output_files(out_dir)
            else:
                # Run locally
                result_files = self.run(fasta_path, out_dir, **kwargs)

            # Parse results into Structure objects
            structures = {}
            for seq_id, pdb_path in result_files.items():
                structure = parser.get_structure(pdb_path)
                structures[seq_id] = structure

            return structures

        finally:
            # Clean up temporary FASTA file
            os.unlink(fasta_path)

    def run(
        self,
        fasta_path: str,
        out_dir: Union[str, Path],
        num_models: int = 5,
        num_recycles: int = 3,
        database_source: Optional[str] = None,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = "alphafold2_ptm",
        db_preset: str = "reduced_dbs",
        additional_options: Optional[List[str]] = None,
    ) -> Dict[str, str]:
        """Execute AlphaFold prediction locally.
        
        Args:
            fasta_path: Path to input FASTA file
            out_dir: Output directory for results
            num_models: Number of models to generate
            num_recycles: Number of recycling iterations
            database_source: Database source for MSA
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
        
        if config.INSTALL_TYPE == "colabfold_container":
            cmd = self._build_colabfold_container_command(
                fasta_path, out_dir, num_models, num_recycles, 
                num_relax, relax_max_iteration, use_templates, 
                model_preset, additional_options
            )
        elif config.INSTALL_TYPE == "alphafold_native":
            cmd = self._build_alphafold_native_command(
                fasta_path, out_dir, max_template_date, 
                model_preset, db_preset, additional_options
            )
        elif config.INSTALL_TYPE == "colabfold_python":
            cmd = self._build_colabfold_python_command(
                fasta_path, out_dir, num_models, num_recycles,
                num_relax, relax_max_iteration, use_templates,
                model_preset, additional_options
            )
        else:
            raise ValueError(f"Unsupported install type: {config.INSTALL_TYPE}")
        
        _LOGGER.info(f"Running AlphaFold command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=out_dir
            )
            _LOGGER.info(f"AlphaFold completed successfully")
            _LOGGER.debug(f"stdout: {result.stdout}")
            
        except subprocess.CalledProcessError as e:
            _LOGGER.error(f"AlphaFold failed with return code {e.returncode}")
            _LOGGER.error(f"stdout: {e.stdout}")
            _LOGGER.error(f"stderr: {e.stderr}")
            raise RuntimeError(f"AlphaFold execution failed: {e.stderr}")
        
        return self._find_output_files(out_dir)

    def make_job(
        self,
        fasta_path: str,
        out_dir: Union[str, Path],
        num_models: int = 5,
        num_recycles: int = 3,
        database_source: Optional[str] = None,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = "alphafold2_ptm",
        db_preset: str = "reduced_dbs",
        additional_options: Optional[List[str]] = None,
        cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None,
    ) -> ClusterJob:
        """Create a cluster job for running AlphaFold.
        
        Args:
            fasta_path: Path to input FASTA file
            out_dir: Output directory for results
            num_models: Number of models to generate
            num_recycles: Number of recycling iterations
            database_source: Database source for MSA
            num_relax: Number of structures to relax
            relax_max_iteration: Maximum relaxation iterations
            use_templates: Whether to use templates
            max_template_date: Maximum template date
            model_preset: Model preset configuration
            db_preset: Database preset
            additional_options: Additional command-line options
            cluster_job_config: Configuration for cluster job
            
        Returns:
            ClusterJob object ready for submission
        """
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)
        
        config = self.config_
        
        if config.INSTALL_TYPE == "colabfold_container":
            cmd = self._build_colabfold_container_command(
                fasta_path, out_dir, num_models, num_recycles, 
                num_relax, relax_max_iteration, use_templates, 
                model_preset, additional_options
            )
        elif config.INSTALL_TYPE == "alphafold_native":
            cmd = self._build_alphafold_native_command(
                fasta_path, out_dir, max_template_date, 
                model_preset, db_preset, additional_options
            )
        elif config.INSTALL_TYPE == "colabfold_python":
            cmd = self._build_colabfold_python_command(
                fasta_path, out_dir, num_models, num_recycles,
                num_relax, relax_max_iteration, use_templates,
                model_preset, additional_options
            )
        else:
            raise ValueError(f"Unsupported install type: {config.INSTALL_TYPE}")
        
        # Create ClusterJob
        if isinstance(cluster_job_config, dict):
            job_config = ClusterJobConfig.from_dict(cluster_job_config)
        else:
            job_config = cluster_job_config or ClusterJobConfig()
        
        # Use ClusterJob.config_job to create the job properly
        if not job_config.has_cluster():
            # Import default cluster
            from enzy_htp.core.clusters import LocalCluster
            job_config.cluster = LocalCluster()
        
        if not job_config.has_res_keywords():
            job_config.res_keywords = {}
        
        job = ClusterJob.config_job(
            commands=' '.join(cmd),
            cluster=job_config.cluster,
            env_settings=[],
            res_keywords=job_config.res_keywords,
            sub_dir=str(out_dir)
        )
        
        return job

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
        
        # Container image
        cmd.append(container_path)
        
        # ColabFold command
        cmd.append("colabfold_batch")
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

    def _build_alphafold_native_command(
        self,
        fasta_path: str,
        out_dir: Path,
        max_template_date: Optional[str],
        model_preset: str,
        db_preset: str,
        additional_options: Optional[List[str]]
    ) -> List[str]:
        """Build command for native AlphaFold execution."""
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

    def _build_colabfold_python_command(
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
        """Build command for ColabFold Python execution."""
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

    def _find_output_files(self, out_dir: Path) -> Dict[str, str]:
        """Find output PDB files in the output directory.
        
        Args:
            out_dir: Output directory to search
            
        Returns:
            Dict mapping sequence IDs to PDB file paths
        """
        result_files = {}
        
        # Look for PDB files
        for pdb_file in out_dir.glob("**/*.pdb"):
            # Extract sequence identifier from filename
            # ColabFold typically names files like: seq_0_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb
            filename = pdb_file.stem
            
            # Try to extract sequence ID (assumes format starts with seq_N)
            if filename.startswith("seq_"):
                seq_id = filename.split("_")[1]
                seq_key = f"seq_{seq_id}"
                
                # Prefer relaxed over unrelaxed, and higher ranked models
                if seq_key not in result_files or (
                    "unrelaxed" in result_files[seq_key] and "unrelaxed" not in str(pdb_file)
                ):
                    result_files[seq_key] = str(pdb_file)
            else:
                # Fallback: use filename as key
                result_files[filename] = str(pdb_file)
        
        return result_files
