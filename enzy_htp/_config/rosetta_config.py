"""Defines RosettaConfig() which holds configuration settings for enzy_htp to interface with the
Rosetta modelling suite. In addition to holding Rosetta. File also contains default_rosetta_config()
which creates a default version of the RosettaConfig() object.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: Qianzhen Shao <shaoqz@icloud.com>
Date: 2023-03-28
"""

import copy
from typing import List, Any, Dict
from copy import deepcopy

from .base_config import BaseConfig
from .armer_config import ARMerConfig

from enzy_htp.core.clusters.accre import Accre


class RosettaConfig(BaseConfig):
    """Class that holds default values for running Rosetta with enzy_htp and 
    also keeps track of required environment variables and exectuables.

    Attributes:
        ROSETTA3 : str() of the environment variable needed for the quite.
        ROSETTA_SCRIPTS : str() of the RosettaScripts executable.
        PARAMS_SCRIPT: str() of the parameter generation file needed for RosettaLigand.
    """

    ROSETTA3: str = "ROSETTA3"
    """Path variable that points to the git installation of rosetta aka 'main'."""

    ROSETTA_SCRIPTS: str = "$ROSETTA3/source/bin/rosetta_scripts.linuxgccrelease"
    """The name of the RosettaScripts executable."""

    PARAMS_SCRIPT: str = f"$ROSETTA3/source/scripts/python/public/molfile_to_params.py"
    """Script used for paramterizing ligands for RosettaLigand protocol."""

    RELAX: str = f"$ROSETTA3/source/bin/relax.default.linuxgccrelease"
    """Executable used to relax a structure/pose."""

    RELAX_MPI: str = f"$ROSETTA3/source/bin/relax.mpi.linuxgccrelease"
    """MPI Executable used to relax a structure/pose."""

    RELAX_MPI_EXEC: str = f"mpiexec -np"
    """MPI Executable used to relax a structure/pose."""

    SCORE: str = f"$ROSETTA3/source/bin/score_jd2.default.linuxgccrelease"
    """Executable used to score a specific structure/pose."""

    PY_2_7: str = "python2.7"
    """Executable or path to the python2.7 executable that Rosetta will use for miscellaneous scripts."""


    DEFAULT_DISTANCE_CONSTRAINT_SETTING:Dict[str, float] = {
        'tolerance': 0.5,
        'penalty': 500.0
    }
    """Default distance constraints used for AtomPair constraints in Rosetta."""

    DEFAULT_ANGLE_CONSTRAINT_SETTING:Dict[str, float] = {
        'tolerance': 20.0,
        'penalty': 50.0
    }
    """Default angle constraints used for Angle constraints in Rosetta."""

    # region == Default values for build_cartesian_ddg_engine() ==
    DEFAULT_CART_DDG_NUM_ITER: int = 10
    """the default value for the NUM_ITER option of cartesian ddg calculation"""

    DEFAULT_CART_DDG_FORCE_ITER: bool = False
    """the default value for the FORCE_ITER option of cartesian ddg calculation"""

    DEFAULT_CART_DDG_SCORE_CUTOFF: float = 1.0
    """the default value for the SCORE_CUTOFF option of cartesian ddg calculation"""

    DEFAULT_CART_DDG_FA_MAX_DIS: float = 9.0
    """the default value for the FA_MAX_DIS option of cartesian ddg calculation"""

    DEFAULT_CART_DDG_SCOREFXN: str = "ref2015_cart"
    """the default value for the SCOREFXN option of cartesian ddg calculation"""

    DEFAULT_CART_DDG_RES_KEYWORDS: Dict = copy.deepcopy(ARMerConfig.SINGLE_CPU_RES) | {'job_name' : 'cart_ddg_EnzyHTP',}
    """The default value for the resource configuration of CART_DDG."""

    DEFAULT_CART_DDG_RELAX_RES_KEYWORDS: Dict = copy.deepcopy(
        ARMerConfig.SINGLE_CPU_RES
    ) | {'job_name' : 'cart_ddg_relax_EnzyHTP',
         'node_cores' : '24',}
    """The default value for the resource configuration of the relax on the WT before CART_DDG."""

    def get_default_cart_ddg_cluster_job_res_keywords(self) -> Dict:
        """function for lazy resolution."""
        return copy.deepcopy(self.DEFAULT_CART_DDG_RES_KEYWORDS)

    def get_default_cart_ddg_cluster_job_config(self) -> Dict:
        """The default value for dictionary that assign arguments to
        ClusterJob.config_job and ClusterJob.wait_to_end during the CART_DDG"""
        return {
            "cluster" : Accre(),
            "res_keywords" : self.get_default_cart_ddg_cluster_job_res_keywords(),
        }

    def get_default_cart_ddg_relax_cluster_job_res_keywords(self) -> Dict:
        """function for lazy resolution."""
        return copy.deepcopy(self.DEFAULT_CART_DDG_RELAX_RES_KEYWORDS)

    def get_default_cart_ddg_relax_cluster_job_config(self) -> Dict:
        """The default value for dictionary that assign arguments to
        ClusterJob.config_job and ClusterJob.wait_to_end during the CART_DDG"""
        return {
            "cluster" : Accre(),
            "res_keywords" : self.get_default_cart_ddg_relax_cluster_job_res_keywords(),
        }

    DEFAULT_CART_DDG_WORK_DIR: str = "./rosetta_cart_ddg"
    """The default value for the cartesian ddg working dir that contains all the temp/result
    files."""
    # endregion

    # region == Default values for relax() ==
    DEFAULT_RELAX_RES_KEYWORDS: Dict = copy.deepcopy(ARMerConfig.SINGLE_CPU_RES) | {
                                            'job_name' : 'rosetta_relax_EnzyHTP',
                                            }
    """The default value for the resource configuration of relax()."""
    
    # endregion

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Rosetta."""
        return [self.ROSETTA_SCRIPTS, self.PARAMS_SCRIPT, self.PY_2_7]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for Rosetta."""
        return [self.ROSETTA3]

    def required_py_modules(self) -> List[str]:
        """ """
        return list()


def default_rosetta_config() -> RosettaConfig:
    """Creates a deep-copied version of the RosettaConfig() class."""
    return deepcopy(RosettaConfig())
