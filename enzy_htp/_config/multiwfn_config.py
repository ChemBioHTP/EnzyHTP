"""Defines a MultiwfnConfig class which is a direct companion to the main MultiwfnInterface class. Stores
the required executables, environment variables and settings for using Multiwfn within enzy_htp. Additionally
features the default_multiwfn_config() method which provides a deep-copied MultiwfnConfig() object with default
values for use in Multiwfn based calculations.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-07-01
"""

from copy import deepcopy
from typing import Any, List, Dict

from .base_config import BaseConfig

from enzy_htp.core.clusters.accre import Accre


class MultiwfnConfig(BaseConfig):
    """Class that holds the default values for running Multiwfn within enzy_htp as well
    the names of required executables and environment variables.

    Attributes:
        EXE : str() corresponding to the Multiwfn executable.
        DIR : str() corresponding to the directory where Multiwfn was installed from.
    """

    EXE: str = "Multiwfn"
    """path of the multiwfn executable"""

    DIR: str = "Multiwfnpath"

    DEFAULT_BOND_DIPOLE_RES_KEYWORDS: Dict = {
        'core_type' : 'cpu',
        'nodes':'1',
        'node_cores' : '4',
        'job_name' : 'bond_dipole_Multiwfn_EnzyHTP',
        'partition' : '<fillthis>',
        'mem_per_core' : '2G', # in GB.
        'walltime' : '5:00:00',
        'account' : '<fillthis>',
        }
    """The default value for the resource configuration of bond dipole calculation."""

    def get_default_bond_dipole_res_keywords(self) -> Dict:
        """function for lazy resolution."""
        return deepcopy(self.DEFAULT_BOND_DIPOLE_RES_KEYWORDS)

    def get_default_bond_dipole_cluster_job_config(self) -> Dict:
        """The default value for dictionary that assign arguments to
        ClusterJob.config_job and ClusterJob.wait_to_end during 
        the bond dipole calculation."""
        return {
            "cluster" : Accre(),
            "res_keywords" : self.get_default_bond_dipole_res_keywords(),
        }

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
