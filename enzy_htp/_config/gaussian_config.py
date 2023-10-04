"""Defines a GaussianConfig class which is a direct companion to the main GaussianInterface class. Stores the
required executables, enviornment variables and settings for using Gaussian within enzy_htp. Additionally
features default_gaussian_config() which provides a deep-copied GaussianConfig() object for use in Gaussian based
calculations.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-11
"""
from copy import deepcopy
from typing import Any, List, Dict

from .base_config import BaseConfig


class GaussianConfig(BaseConfig):
    """Class that holds default values for running Gaussian within enzy_htp as well as the required
    executables and variables that need to be present in the environment.


    Attributes:
        G16_EXE : str() corresponding to the Gaussian 16 executable.
        G09_EXE	: str() corresponding to the Gaussian 09 executable.
        OM_LVL : list() holding default values for ONIOM simulation header.
        ONIOM_KW : str() with line of default ONIOM keywords.
        KEYWORDS : dict() holding settings for ONIOM simulation.
        LAYER_PRESET : int() default ONIOM level for simulation.
        LAYER_ATOMS : list() holding the atom names.
    """

    G16_EXE: str = "g16"
    """Variable for the g16 version of Gaussian."""

    G09_EXE: str = "g09"
    """Variable for the g09 version of Gaussian."""

    OM_LVL: List[str] = [
        "wb97xd/6-31g(d)",
        "amber=hardfirst",
    ]
    """Keywords for an ONIOM job. Use mechanic embedding as default / hardfirst means find user assigned parameters first / iops are for reduce the size of output."""

    ONIOM_KW: str = f"oniom({':'.join(OM_LVL)})"
    """Default ONIOM keywords. Note that it DOES NOT update when the OM_LEVEL keywords are changed."""

    KEYWORDS: Dict[str, List[str]] = {
        "spe": [
            ONIOM_KW,
            "nosymm",
            "geom=connectivity",
            "iop(2/9=1000,3/33=0,3/36=-1,3/160=2,4/33=0,5/33=0,6/7=2,6/9=2,6/10=2,6/12=2,7/8=10,2/15=3)",
        ],
        "opt": [
            ONIOM_KW,
            "opt=(modredundant,calcfc)",
            "freq=(selectnormalmodes,hpmodes)",
            "nosymm",
            "geom=connectivity",
            "iop(2/9=1000,3/33=0,3/36=-1,3/160=2,4/33=0,5/33=0,6/7=2,6/9=2,6/10=2,6/12=2,7/8=10,2/15=3)",
        ],
        "tsopt": [
            ONIOM_KW,
            "opt=(calcfc,ts,noeigen,gdiis)",
            "freq=(selectnormalmodes,hpmodes)",
            "nosymm",
            "geom=connectivity",
            "iop(2/9=1000,3/33=0,3/36=-1,3/160=2,4/33=0,5/33=0,6/7=2,6/9=2,6/10=2,6/12=2,7/8=10,2/15=3)",
        ],
    }

    LAYER_PRESET: int = 0
    """Default layer preset for the ONIOM simulation."""

    LAYER_ATOMS: List[str] = []
    """Default Atoms in the ONIOM layer."""

    def __init__(self, parent=None):
        self._parent = parent

    def required_executables(self) -> List[str]:
        """ """
        return [self.G16_EXE]

    def required_env_vars(self):
        """ """
        return []

    def required_py_modules(self) -> List[str]:
        """ """
        return list()

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


def default_gaussian_config() -> GaussianConfig:
    """Creates a deep-copied default version of the GaussianConfig() class."""
    return deepcopy(GaussianConfig())
