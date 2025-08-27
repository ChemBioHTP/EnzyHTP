"""Defines AlphafoldConfig() which holds configuration settings for enzy_htp to interface with the 
AlphaFold software package.

Author: Gemini
Date: 2025-08-16
"""
from __future__ import annotations
import os
from typing import List, Dict
from copy import deepcopy
from enzy_htp.core.general import get_str_for_print_class_var
from enzy_htp.core.logger import _LOGGER

from .base_config import BaseConfig


class AlphafoldConfig(BaseConfig):
    """Class that holds default values for running AlphaFold2 within enzy_htp."""
    
    # Installation configuration
    INSTALL_TYPE: str = "colabfold_container"
    """Installation type for AlphaFold2.
    Options: 'colabfold_container', 'alphafold2_native_container', 'alphafold2_native_python'
    """
    
    CONTAINER_TYPE: str = "apptainer"
    """Container software type.
    Options: 'docker', 'apptainer', 'singularity'
    """
    
    CONTAINER_PATH: str = "~/bin/colabfold_1.5.5-cuda12.2.2.sif"
    """Path to the container image file."""
    
    CONTAINER_BIND_PATHS: dict = {
        "/cache": "~/bin/colabfold/cache"
    }
    """Bind paths for container mounting.
    Maps container paths to host paths. Work directory is handled dynamically.
    """
    
    # Native AlphaFold2 configuration
    EXECUTABLE_PATH: str = "colabfold_batch"
    """Path to AlphaFold2 executable.
    For ColabFold: 'colabfold_batch'
    For native AlphaFold2: path to run_alphafold.py or run_docker.py
    """
    
    DATA_DIR: str = ""
    """Path to AlphaFold2 data directory.
    Required for native AlphaFold2 installations.
    """
    
    # Database paths for native AlphaFold2 
    UNIREF90_DATABASE_PATH: str = ""
    """Path to UniRef90 database."""
    
    MGNIFY_DATABASE_PATH: str = ""
    """Path to MGnify database."""
    
    UNIREF30_DATABASE_PATH: str = ""
    """Path to UniRef30 database."""
    
    BFD_DATABASE_PATH: str = ""
    """Path to BFD database."""
    
    PDB70_DATABASE_PATH: str = ""
    """Path to PDB70 database."""
    
    TEMPLATE_MMCIF_DIR: str = ""
    """Path to template mmCIF directory."""
    
    PDB_SEQRES_DATABASE_PATH: str = ""
    """Path to PDB seqres database."""
    
    OBSOLETE_PDBS_PATH: str = ""
    """Path to obsolete PDBs file."""
    
    UNIPROT_DATABASE_PATH: str = ""
    """Path to UniProt database."""
    
    USE_GPU_RELAX: bool = True
    """Whether to use GPU for relaxation."""
    
    # Common configuration
    CACHE_DIR: str = "~/bin/colabfold/cache"
    """Path to cache directory for MSA and model downloads."""
    
    WORK_DIR: str = "./alphafold2_predictions"
    """Default working directory for predictions."""

    # Default resource settings for cluster jobs
    DEFAULT_AF2_CLUSTER_JOB_RES_KEYWORDS = {
        "gpu": {
            'core_type': 'gpu',
            'nodes': '1',
            'node_cores': 'nvidia_rtx_a4000:1',
            'job_name': 'AF2_EnzyHTP',
            'partition': '<fillthis>',
            'account': '<fillthis>',
            'mem_per_core': '24G',
            'walltime': '16:00:00',
        },
        "cpu": {
            'core_type': 'cpu',
            'nodes': '1',
            'node_cores': '24',
            'job_name': 'AF2_EnzyHTP',
            'partition': '<fillthis>',
            'account': '<fillthis>',
            'mem_per_core': '2G',
            'walltime': '16:00:00',
        }
    }
    """Default res_keywords for AlphaFold2 cluster jobs."""

    def get_default_af2_cluster_job_res_keywords(self, key: str) -> Dict:
        """Get default resource keywords for AF2 cluster jobs."""
        return deepcopy(self.DEFAULT_AF2_CLUSTER_JOB_RES_KEYWORDS[key])

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Alphafold."""
        if self.INSTALL_TYPE == "colabfold_container":
            return [self.CONTAINER_TYPE]
        return []

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for Alphafold."""
        return []

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required python modules for Alphafold."""
        return []

    @classmethod
    def display(cls) -> None:
        """Method that prints out all settings for the AlphafoldConfig class to the stdout."""
        dash_line: str = "-" * 40
        dis_info = f"AlphafoldConfig settings:{os.linesep}{dash_line}{get_str_for_print_class_var(cls)}"
        _LOGGER.info(dis_info)

def default_alphafold_config() -> AlphafoldConfig:
    """Creates a deep-copied default version of the AlphafoldConfig() class."""
    return AlphafoldConfig()