"""TODO

Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2023-04-02
"""
from typing import List

from copy import deepcopy

from .base_config import BaseConfig


class MOEConfig(BaseConfig):
    """Class that holds default values and environment requirements for running MOE within
    enzy_htp.

    Attributes:
        MOE : str() environment variable referencing the home directory of the package.
        MOE_BATCH : str() of the name of the moebatch executable.

    """

    MOE: str = "MOE"
    """MOE environment variable which points to the """

    MOE_BATCH: str = "$MOE/bin/moebatch"
    """Path of the moebatch executable which can run MOE functionality at the commandline."""

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for the MOE."""
        return [self.MOE_BATCH]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for the MOE."""
        return [self.MOE]

    def required_py_modules(self) -> List[str]:
        """ """
        return list()


def default_moe_config() -> MOEConfig:
    """Creates a deep-copied default version of the MOEConfig() class."""
    return deepcopy(MOEConfig())
