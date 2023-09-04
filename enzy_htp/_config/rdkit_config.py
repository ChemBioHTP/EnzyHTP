""" """
#TODO(CJ): do the documentation here
from copy import deepcopy
from typing import List

from .base_config import BaseConfig

class RDKitConfig(BaseConfig):

    SUPPORTED_FTYPES:List[str] = ".pdb .mol .mol2 .sdf".split()
    """ """

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for xtb."""
        return list()

    def required_env_vars(self)-> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return list()


    def required_py_modules(self)-> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return ["rdkit"]

def default_rdkit_config() -> RDKitConfig:
    """Creates a deep-copied default version of the RDKitConfig() class."""
    return deepcopy(RDKitConfig())
