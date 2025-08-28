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
    
    DATA_DIR: str = "/sb/apps/alphafold-data.230"
    """Path to AlphaFold2 data directory.
    Required for native AlphaFold2 installations.
    """
    
    # Database paths for native AlphaFold2 
    UNIREF90_DATABASE_PATH: str = "{DATA_DIR}/uniref90/uniref90.fasta"
    """Path to UniRef90 database."""
    
    MGNIFY_DATABASE_PATH: str = "{DATA_DIR}/mgnify/mgy_clusters_2022_05.fa"
    """Path to MGnify database."""
    
    UNIREF30_DATABASE_PATH: str = "{DATA_DIR}/uniref30/UniRef30_2021_03"
    """Path to UniRef30 database."""
    
    BFD_DATABASE_PATH: str = "{DATA_DIR}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    """Path to BFD database."""
    
    PDB70_DATABASE_PATH: str = "{DATA_DIR}/pdb70/pdb70"
    """Path to PDB70 database."""
    
    TEMPLATE_MMCIF_DIR: str = "{DATA_DIR}/pdb_mmcif/mmcif_files"
    """Path to template mmCIF directory."""
    
    PDB_SEQRES_DATABASE_PATH: str = "{DATA_DIR}/pdb_seqres/pdb_seqres.txt"
    """Path to PDB seqres database."""
    
    OBSOLETE_PDBS_PATH: str = "{DATA_DIR}/pdb_mmcif/obsolete.dat"
    """Path to obsolete PDBs file."""
    
    UNIPROT_DATABASE_PATH: str = "{DATA_DIR}/uniprot/uniprot.fasta"
    """Path to UniProt database."""
    
    SMALL_BFD_DATABASE_PATH: str = "{DATA_DIR}/small_bfd/bfd-first_non_consensus_sequences.fasta"
    """Path to small BFD database for reduced_dbs preset."""
    
    # Binary paths for native AlphaFold2
    ALPHAFOLD_BIN_DIR: str = "/sb/apps/alphafold232/miniconda3/envs/af232/bin"
    """Base directory for AlphaFold2 binary executables."""
    
    HHBLITS_BINARY_PATH: str = "{ALPHAFOLD_BIN_DIR}/hhblits"
    """Path to HHblits executable."""
    
    HHSEARCH_BINARY_PATH: str = "{ALPHAFOLD_BIN_DIR}/hhsearch"
    """Path to HHsearch executable."""
    
    HMMBUILD_BINARY_PATH: str = "{ALPHAFOLD_BIN_DIR}/hmmbuild"
    """Path to hmmbuild executable."""
    
    HMMSEARCH_BINARY_PATH: str = "{ALPHAFOLD_BIN_DIR}/hmmsearch"
    """Path to hmmsearch executable."""
    
    JACKHMMER_BINARY_PATH: str = "{ALPHAFOLD_BIN_DIR}/jackhmmer"
    """Path to JackHMMER executable."""
    
    KALIGN_BINARY_PATH: str = "{ALPHAFOLD_BIN_DIR}/kalign"
    """Path to Kalign executable."""
    
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
    
    # Database path getter methods
    def get_uniref90_database_path(self) -> str:
        """Get the UniRef90 database path, expanding DATA_DIR if needed."""
        if self.UNIREF90_DATABASE_PATH and "{DATA_DIR}" in self.UNIREF90_DATABASE_PATH:
            return self.UNIREF90_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.UNIREF90_DATABASE_PATH
    
    def get_mgnify_database_path(self) -> str:
        """Get the MGnify database path, expanding DATA_DIR if needed."""
        if self.MGNIFY_DATABASE_PATH and "{DATA_DIR}" in self.MGNIFY_DATABASE_PATH:
            return self.MGNIFY_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.MGNIFY_DATABASE_PATH
    
    def get_uniref30_database_path(self) -> str:
        """Get the UniRef30 database path, expanding DATA_DIR if needed."""
        if self.UNIREF30_DATABASE_PATH and "{DATA_DIR}" in self.UNIREF30_DATABASE_PATH:
            return self.UNIREF30_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.UNIREF30_DATABASE_PATH
    
    def get_bfd_database_path(self) -> str:
        """Get the BFD database path, expanding DATA_DIR if needed."""
        if self.BFD_DATABASE_PATH and "{DATA_DIR}" in self.BFD_DATABASE_PATH:
            return self.BFD_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.BFD_DATABASE_PATH
    
    def get_pdb70_database_path(self) -> str:
        """Get the PDB70 database path, expanding DATA_DIR if needed."""
        if self.PDB70_DATABASE_PATH and "{DATA_DIR}" in self.PDB70_DATABASE_PATH:
            return self.PDB70_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.PDB70_DATABASE_PATH
    
    def get_template_mmcif_dir(self) -> str:
        """Get the template mmCIF directory path, expanding DATA_DIR if needed."""
        if self.TEMPLATE_MMCIF_DIR and "{DATA_DIR}" in self.TEMPLATE_MMCIF_DIR:
            return self.TEMPLATE_MMCIF_DIR.format(DATA_DIR=self.DATA_DIR)
        return self.TEMPLATE_MMCIF_DIR
    
    def get_pdb_seqres_database_path(self) -> str:
        """Get the PDB seqres database path, expanding DATA_DIR if needed."""
        if self.PDB_SEQRES_DATABASE_PATH and "{DATA_DIR}" in self.PDB_SEQRES_DATABASE_PATH:
            return self.PDB_SEQRES_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.PDB_SEQRES_DATABASE_PATH
    
    def get_obsolete_pdbs_path(self) -> str:
        """Get the obsolete PDBs file path, expanding DATA_DIR if needed."""
        if self.OBSOLETE_PDBS_PATH and "{DATA_DIR}" in self.OBSOLETE_PDBS_PATH:
            return self.OBSOLETE_PDBS_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.OBSOLETE_PDBS_PATH
    
    def get_uniprot_database_path(self) -> str:
        """Get the UniProt database path, expanding DATA_DIR if needed."""
        if self.UNIPROT_DATABASE_PATH and "{DATA_DIR}" in self.UNIPROT_DATABASE_PATH:
            return self.UNIPROT_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.UNIPROT_DATABASE_PATH
    
    def get_small_bfd_database_path(self) -> str:
        """Get the small BFD database path, expanding DATA_DIR if needed."""
        if self.SMALL_BFD_DATABASE_PATH and "{DATA_DIR}" in self.SMALL_BFD_DATABASE_PATH:
            return self.SMALL_BFD_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)
        return self.SMALL_BFD_DATABASE_PATH
    
    # Binary path getter methods
    def get_hhblits_binary_path(self) -> str:
        """Get the HHblits binary path, expanding ALPHAFOLD_BIN_DIR if needed."""
        if self.HHBLITS_BINARY_PATH and "{ALPHAFOLD_BIN_DIR}" in self.HHBLITS_BINARY_PATH:
            return self.HHBLITS_BINARY_PATH.format(ALPHAFOLD_BIN_DIR=self.ALPHAFOLD_BIN_DIR)
        return self.HHBLITS_BINARY_PATH
    
    def get_hhsearch_binary_path(self) -> str:
        """Get the HHsearch binary path, expanding ALPHAFOLD_BIN_DIR if needed."""
        if self.HHSEARCH_BINARY_PATH and "{ALPHAFOLD_BIN_DIR}" in self.HHSEARCH_BINARY_PATH:
            return self.HHSEARCH_BINARY_PATH.format(ALPHAFOLD_BIN_DIR=self.ALPHAFOLD_BIN_DIR)
        return self.HHSEARCH_BINARY_PATH
    
    def get_hmmbuild_binary_path(self) -> str:
        """Get the hmmbuild binary path, expanding ALPHAFOLD_BIN_DIR if needed."""
        if self.HMMBUILD_BINARY_PATH and "{ALPHAFOLD_BIN_DIR}" in self.HMMBUILD_BINARY_PATH:
            return self.HMMBUILD_BINARY_PATH.format(ALPHAFOLD_BIN_DIR=self.ALPHAFOLD_BIN_DIR)
        return self.HMMBUILD_BINARY_PATH
    
    def get_hmmsearch_binary_path(self) -> str:
        """Get the hmmsearch binary path, expanding ALPHAFOLD_BIN_DIR if needed."""
        if self.HMMSEARCH_BINARY_PATH and "{ALPHAFOLD_BIN_DIR}" in self.HMMSEARCH_BINARY_PATH:
            return self.HMMSEARCH_BINARY_PATH.format(ALPHAFOLD_BIN_DIR=self.ALPHAFOLD_BIN_DIR)
        return self.HMMSEARCH_BINARY_PATH
    
    def get_jackhmmer_binary_path(self) -> str:
        """Get the JackHMMER binary path, expanding ALPHAFOLD_BIN_DIR if needed."""
        if self.JACKHMMER_BINARY_PATH and "{ALPHAFOLD_BIN_DIR}" in self.JACKHMMER_BINARY_PATH:
            return self.JACKHMMER_BINARY_PATH.format(ALPHAFOLD_BIN_DIR=self.ALPHAFOLD_BIN_DIR)
        return self.JACKHMMER_BINARY_PATH
    
    def get_kalign_binary_path(self) -> str:
        """Get the Kalign binary path, expanding ALPHAFOLD_BIN_DIR if needed."""
        if self.KALIGN_BINARY_PATH and "{ALPHAFOLD_BIN_DIR}" in self.KALIGN_BINARY_PATH:
            return self.KALIGN_BINARY_PATH.format(ALPHAFOLD_BIN_DIR=self.ALPHAFOLD_BIN_DIR)
        return self.KALIGN_BINARY_PATH

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Alphafold."""
        if self.INSTALL_TYPE == "colabfold_container":
            return [self.CONTAINER_TYPE]
        if self.INSTALL_TYPE == "alphafold2_native_python":
            # Return executables that are actually configured using getter functions
            executables = []
            hhblits_path = self.get_hhblits_binary_path()
            if hhblits_path:
                executables.append(hhblits_path)
            hhsearch_path = self.get_hhsearch_binary_path()
            if hhsearch_path:
                executables.append(hhsearch_path)
            hmmbuild_path = self.get_hmmbuild_binary_path()
            if hmmbuild_path:
                executables.append(hmmbuild_path)
            hmmsearch_path = self.get_hmmsearch_binary_path()
            if hmmsearch_path:
                executables.append(hmmsearch_path)
            jackhmmer_path = self.get_jackhmmer_binary_path()
            if jackhmmer_path:
                executables.append(jackhmmer_path)
            kalign_path = self.get_kalign_binary_path()
            if kalign_path:
                executables.append(kalign_path)
            return executables
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