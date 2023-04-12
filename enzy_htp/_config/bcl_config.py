"""TODO

Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2023-04-02
"""
from typing import List

from copy import deepcopy


class BCLConfig:

    BCL_EXE:str = "bcl.exe"
    """The main bcl executable."""

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for the BCL."""
        return [
            self.BCL_EXE
        ]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for the BCL."""
        return []




def default_bcl_config() -> BCLConfig:
    """Creates a deep-copied default version of the BCLConfig() class."""
    return deepcopy(BCLConfig())
