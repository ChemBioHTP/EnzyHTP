"""Defines a GaussianConfig class which is a direct companion to the main GaussianInterface class. Stores the
required executables, enviornment variables and settings for using Gaussian within enzy_htp. Additionally
features default_gaussian_config() which provides a deep-copied GaussianConfig() object for use in Gaussian based
calculations.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-11
"""
from copy import deepcopy
import copy
from typing import Any, List, Dict

from .base_config import BaseConfig
from .armer_config import ARMerConfig

from enzy_htp.chemical import QMLevelOfTheory, MMLevelOfTheory
from enzy_htp.core.clusters.accre import Accre

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
    DEFAULT_SOLVENT_MODEL: str = "SMD"
    """default solvent model when a solvent is specificed but not solvent model"""

    GAUSSIAN_EXE: str = "g16"
    """Variable for the g16 version of Gaussian."""

    # region == Default values for build_single_point_engine() ==
    DEFAULT_SPE_NAME: str = "gaussian_spe"
    """The default value for the name tag of the SPE engine"""

    DEFAULT_SPE_METHOD: QMLevelOfTheory = QMLevelOfTheory(
        basis_set = "def2-svp",
        method = "pbe0",
        solvent = None,
        solv_method = None,
    )
    """The default value for the level of theory of the SPE engine"""

    DEFAULT_SPE_KEEP_GEOM: bool = True
    """The default value for whether the SPE engine keep/restore the geometry
    orientation. This triggers the use of nosymm for now. TODO probably better
    to restore instead."""

    DEFAULT_SPE_CAPPING_METHOD: str = "h_cap"
    """The default value for the capping method of the SPE calculation
    when a subregion is specified."""

    DEFAULT_SPE_RES_KEYWORDS: Dict = copy.deepcopy(ARMerConfig.QM_SPE_CPU_RES)
    """The default value for the resource configuration of SPE."""

    def get_default_qm_spe_cluster_job_res_keywords(self) -> Dict:
        """function for lazy resolution."""
        return copy.deepcopy(self.DEFAULT_SPE_RES_KEYWORDS)

    def get_default_qm_spe_cluster_job(self) -> Dict:
        """The default value for dictionary that assign arguments to
        ClusterJob.config_job and ClusterJob.wait_to_end during the QM SPE"""
        return {
            "cluster" : Accre(),
            "res_keywords" : self.get_default_qm_spe_cluster_job_res_keywords(),
        }

    DEFAULT_SPE_WORK_DIR: str = "./QM_SPE"
    """The default value for the SPE working dir that contains all the temp/result
    files."""
    # endregion

    # region == TODO ==
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
    # endregion == TODO ==

    def __init__(self, parent=None):
        self._parent = parent

    def required_executables(self) -> List[str]:
        """ """
        return [self.GAUSSIAN_EXE, "formchk", "cubegen"]

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
