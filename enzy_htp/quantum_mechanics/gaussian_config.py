"""Defines a GaussianConfig class that serves as a bridge for enzy_htp to utilize Gaussian software.

TODO(CJ): elaborate this part
Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-11
"""

from copy import deepcopy
from typing import Any, List, Dict


class GaussianConfig:
    """Class that holds default values for running Gaussian within enzy_htp and al

    Attributes:
	    G16_EXE : str() corresponding to the Gaussian 16 executable.	
	    G09_EXE	: str() corresponding to the Gaussian 09 executable.
		OM_LVL : list() holding default values for ONIOM simulation header.
		ONIOM_KW : str() with line of default ONIOM keywords.
		KEYWORDS : dict() holding settings for ONIOM simulation.
        LAYER_PRESET : int() default ONIOM level for simulation.
		LAYER_ATOMS : list() holding the atom names.
    """

    G16_EXE:str = "g16"
    """Variable for the g16 version of Gaussian."""

    G09_EXE:str = "g09"
    """Variable for the g09 version of Gaussian."""

    OM_LVL:List[str] = [
        "wb97xd/6-31g(d)",
        "amber=hardfirst",
    ]
    """Keywords for an ONIOM job. Use mechanic embedding as default / hardfirst means find user assigned parameters first / iops are for reduce the size of output."""

    ONIOM_KW:str = f"oniom({':'.join(OM_LVL)})"
    """Default ONIOM keywords. Note that it DOES NOT update when the OM_LEVEL keywords are changed."""

    KEYWORDS:Dict[str,List[str]] = {
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

    LAYER_PRESET:int = 0
    """Default layer preset for the ONIOM simulation."""
    
    LAYER_ATOMS:List[str] = []
    """Default Atoms in the ONION layer."""

    def __init__(self, parent=None):
        self._parent = parent

    def required_executables(self) -> List[str]:
        return [self.G16_EXE, self.G09_EXE]

    def required_env_vars(self):
        return []

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        is_path = False
        if key in {"G09_EXE", "G16_EXE"}:
            is_path = True

        setattr(self, key, value)
        if is_path:
            self._parent.update_paths()


def default_gaussian_config() -> GaussianConfig:
    """Creates a deep-copied default version of the GaussianConfig() class."""
    return deepcopy(GaussianConfig())
