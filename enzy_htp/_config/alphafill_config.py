"""Defines an AlphaFillConfig() which holds configuration settings for enzy_htp to interface with the AlphaFill 
software package. File also contains default_alphafill_config() which creates a default version of 
the AlphaFillConfig() class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-09-14
"""
from typing import List, Any

from copy import deepcopy

from .base_config import BaseConfig


class AlphaFillConfig(BaseConfig):
    """Class that holds default values for running AlphaFill within enzy_htp.

    Attributes:
        ALPHAFILL_EXE: The path to the AlphaFill .exe as a str().
        CONFIG_FILE: The path to the AlphaFill .conf file as a str().
    """

    ALPHAFILL_EXE: str = "alphafill"
    """Path to the AlphaFill .exe as a str()."""

    CONFIG_FILE: str = "alphafill.conf"
    """Path to the AlphaFill .conf file as a str()."""

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for AlphaFill."""
        return [self.ALPHAFILL_EXE]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for AlphaFill."""
        return []

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required python modules for AlphaFill."""
        return list()


def default_alphafill_config() -> AlphaFillConfig:
    """Creates a deep-copied default version of the AlphaFillConfig() class."""
    return deepcopy(AlphaFillConfig())
