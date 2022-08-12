"""Defines a MultiwfnConfig class which is a direct companion to the main MultiwfnInterface class. Stores
the required executables, environment variables and settings for using Multiwfn within enzy_htp. Additionally
features the default_multiwfn_config() method which provides a deep-copied MultiwfnConfig() object with default
values for use in Multiwfn based calculations.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-07-01
"""

from copy import deepcopy
from typing import Any, List


class MultiwfnConfig:
    """Class that holds the default values for running Multiwfn within enzy_htp as well
    the names of required executables and environment variables.

    Attributes:
        EXE : str() corresponding to the Multiwfn executable.
        DIR : str() corresponding to the directory where Multiwfn was installed from.
    """

    EXE: str = "Multiwfn"
    """Name of the """

    DIR: str = "$Multiwfnpath"

    def required_executables(self):
        return [self.EXE]

    def required_env_vars(self):
        return [self.DIR]

    def __getitem__(self, key: str) -> Any:
        if key.count("."):
            key1, key2 = key.split(".")
            return getattr(self, key1)[key2]
        else:
            return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        if key.count("."):
            key1, key2 = key.split(".")
            GaussianConfig.__dict__[key1][key2] = value
        else:
            setattr(self, key, value)


def default_multiwfn_config() -> MultiwfnConfig:
    """Creates a deep-copied default version of the MultiwfnConfig() class."""
    return deepcopy(MultiwfnConfig())
