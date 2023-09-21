"""Defines a MultiwfnConfig class which is a direct companion to the main MultiwfnInterface class. Stores
the required executables, environment variables and settings for using Multiwfn within enzy_htp. Additionally
features the default_multiwfn_config() method which provides a deep-copied MultiwfnConfig() object with default
values for use in Multiwfn based calculations.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-07-01
"""

from copy import deepcopy
from typing import Any, List

from .base_config import BaseConfig


class MultiwfnConfig(BaseConfig):
    """Class that holds the default values for running Multiwfn within enzy_htp as well
    the names of required executables and environment variables.

    Attributes:
        EXE : str() corresponding to the Multiwfn executable.
        DIR : str() corresponding to the directory where Multiwfn was installed from.
    """

    EXE: str = "Multiwfn"
    """Name of the """

    DIR: str = "Multiwfnpath"

    def required_executables(self) -> List[str]:
        """ """
        return [self.EXE]

    def required_env_vars(self) -> List[str]:
        """ """
        return [self.DIR]

    def required_py_modules(self) -> List[str]:
        """ """
        return list()


def default_multiwfn_config() -> MultiwfnConfig:
    """Creates a deep-copied default version of the MultiwfnConfig() class."""
    return deepcopy(MultiwfnConfig())
