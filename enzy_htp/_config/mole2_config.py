"""Defines an Mole2Config() which holds configuration settings for enzy_htp to interface with the Mole2 
software package. File also contains default_mole2_config() which creates a default version of 
the Mole2Config() class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-09-26
"""

from typing import List, Any

from copy import deepcopy

from .base_config import BaseConfig


class Mole2Config(BaseConfig):
    """Class that holds default values for running Mole2 within enzy_htp.

    Attributes:
        PROBE: The probe radius to use in a run in A.
        INNER: The inner radius to use in a run in A.
        MESH_DENSITY: The mesh density to use in a run in A.
        MONO: The path to mono as a str().
        MOLE2: The path to Mole2 as a str().
        IGNORE_HETATM: Should HETATM records be considered?
    """

    IGNORE_HETATM:bool=False 
    """Should HETATM records be considered?"""

    PROBE:float=3.5
    """The probe radius to use in a run in A."""

    INNER:float=1.5
    """The inner radius to use in a run in A."""

    MESH_DENSITY:float=0.5
    """The mesh density to use in a run in A."""

    MONO:str="mono"
    """Path to mono as a str()."""

    MOLE2:str="~/Downloads/Mole2_plugin/MOLE25_binary/mole2.exe"
    """Path to Mole2 .exe as a str()."""


    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Mole2."""
        return [self.MOLE2, self.MONO ]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for Mole2."""
        return []

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required python modules for Mole2."""
        return list()


def default_mole2_config() -> Mole2Config:
    """Creates a deep-copied default version of the Mole2Config() class."""
    return deepcopy(Mole2Config())
