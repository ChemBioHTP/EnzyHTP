"""Defines RDKitConfig() which holds confirguration settings for enzy_htp to interface with the 
RDKit software package. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2025-06-27
"""
from copy import deepcopy
from typing import List

from .base_config import BaseConfig


class RDKitConfig(BaseConfig):
    """Class that holds configuration settings for RDKit, including supported file types and
    the required executables, environment variables, and Python modules.

    Attributes:
        SUPPORTED_FTYPES : List[str] holding supported file types for RDKit operations.
    """

    SUPPORTED_FTYPES: List[str] = ".pdb .mol .mol2 .sdf".split()
    """Supported file types for RDKit operations."""

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for xtb."""
        return list()

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return list()

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return ["rdkit"]


def default_rdkit_config() -> RDKitConfig:
    """Creates a deep-copied default version of the RDKitConfig() class."""
    return deepcopy(RDKitConfig())
