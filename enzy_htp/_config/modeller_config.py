"""*NOTE* THIS MODULE IS DEPROCATED FOR NOW

Defines ModellerConfig() which holds configuration settings for enzy_htp to interface with the
Modeller modelling suite. In addition to holding Modeller. File also contains default_rosetta_config()
which creates a default version of the ModellerConfig() object.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-05-30
"""
from typing import List, Any

from copy import deepcopy

from .base_config import BaseConfig


class ModellerConfig(BaseConfig):
    """Class that holds default values for running the Modeller python package."""

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for AlphaFill."""
        return []

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for AlphaFill."""
        return []

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required python modules for AlphaFill."""
        return ["modeller", "modeller.automodel"]

def default_modeller_config() -> ModellerConfig:
    """Creates a deep-copied default version of the ModellerConfig() class."""

    return deepcopy(ModellerConfig())
