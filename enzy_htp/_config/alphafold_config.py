"""Defines AlphafoldConfig() which holds configuration settings for enzy_htp to interface with the 
AlphaFold software package.

Author: Gemini
Date: 2025-08-16
"""
from __future__ import annotations
import os
from typing import List
from enzy_htp.core.general import get_str_for_print_class_var
from enzy_htp.core.logger import _LOGGER

from .base_config import BaseConfig


class AlphafoldConfig(BaseConfig):
    """Class that holds default values for running Alphafold within enzy_htp.

    Attributes:
        INSTALL_TYPE: str, installation type ('colabfold_container', 'alphafold_native_container', 'alphafold_native_python')
        CONTAINER_TYPE: str, container software type ('docker', 'apptainer', 'singularity')
        CONTAINER_PATH: str, path to the container image
        CONTAINER_BIND_PATHS: dict, bind paths for container
        EXECUTABLE_PATH: str, path to alphafold executable
        DATA_DIR: str, path to alphafold data directory
        CACHE_DIR: str, path to cache directory
        WORK_DIR: str, default working directory for predictions
    """
    
    # Installation configuration
    INSTALL_TYPE: str = "colabfold_container"  # colabfold_container, alphafold_native_container, alphafold_native_python
    CONTAINER_TYPE: str = "apptainer"  # docker, apptainer, singularity
    CONTAINER_PATH: str = "~/bin/colabfold_1.5.5-cuda12.2.2.sif"
    CONTAINER_BIND_PATHS: dict = {
        "/cache": "~/bin/colabfold/cache",
        "/work": "~/colabfold-test/"
    }
    
    # Native AlphaFold configuration
    EXECUTABLE_PATH: str = "colabfold_batch"  # For colabfold, or path to run_alphafold.py for native
    DATA_DIR: str = ""  # Required for native alphafold
    
    # Common configuration
    CACHE_DIR: str = "~/bin/colabfold/cache"
    WORK_DIR: str = "./alphafold_predictions"

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