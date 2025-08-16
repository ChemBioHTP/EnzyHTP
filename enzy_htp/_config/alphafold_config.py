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
        HOME : str() corresponding to Alphafold home directory on the system.
    """

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Alphafold."""
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