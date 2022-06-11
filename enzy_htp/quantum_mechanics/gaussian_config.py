"""Defines a GaussianConfig class that serves as a bridge for enzy_htp to utilize Gaussian software.

TODO(CJ): elaborate this part
Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-11
"""

from copy import deepcopy
from typing import Any, List


class GaussianConfig:
    # -----------------------------
    #   >>>>>>>>ONIOM<<<<<<<<
    # -----------------------------
    def __init__(self, parent=None):
        self._parent = parent

    def required_executables(self) -> List[str]:
        return [self.G16_EXE, self.G09_EXE]

    def required_env_vars(self):
        return []

    # -----------------------------
    # Cores for gaussian job (higher pirority)
    #
    # n_cores = n_cores
    # -----------------------------
    # Per core memory in MB for gaussian job (higher pirority)
    #
    # max_core = max_core
    # -----------------------------
    # Executable g16 command for current environment
    #
    G16_EXE = "g16"
    # -----------------------------
    # Executable g16 command for current environment
    #
    G09_EXE = "g09"
    # -----------------------------
    # key words for ONIOM job (strategy: use work type as a key to indicate different keyword set.)
    # Use mechanic embedding as default / hardfirst means find user assigned parameters first / iops are for reduce the size of output.
    #
    # for oniom method level
    # -- open for edit --
    OM_LVL = [
        "wb97xd/6-31g(d)",
        "amber=hardfirst",
    ]  # add in the order of theory level  e.g.: ['wb97xd/def2svp', 'PM7', 'amber']
    # -- open for edit --
    ONIOM_KW = "oniom(" + ":".join(OM_LVL) + ")"

    # complete keywords
    KEYWORDS = {
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
    # -----------------------------
    # layer settings for oniom
    # preset: 0: no preset (fill layer atom manually) 1: preset_1 -> xxx
    LAYER_PRESET = 0
    LAYER_ATOMS = []

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
