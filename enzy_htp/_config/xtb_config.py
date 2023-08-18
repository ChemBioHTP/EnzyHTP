"""Defines a  XTBConfig() which holds configuration settings for enzy_htp to interface with the xtb
software package. File also contains default_bcl_config() which creates a default version of 
the XTBConfig() class.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-08-17
"""
#TODO(CJ): update this stuff
from typing import List
from copy import deepcopy

from .base_config import BaseConfig

class XTBConfig(BaseConfig):
    """Class that holds default alues for running xtb within enzy_htp.

    Attributes:
        XTB_EXE : str() corresponding to xtb application.
        N_ITER : int() corresponding to number of scf iterations to run.
    """

    XTB_EXE:str = "xtb"
    """Name of the xtb executable as a str()."""

    N_ITER:int = 500
    """Number of SCF iterations to try as an int()."""

    N_PROC:int = 1 
    """Number of parrallel processes to use during a calculation as an int()."""

    KMP_STACKSIZE:str = "KMP_STACKSIZE"
    """Environment variable corresponding to KMP stack size."""

    OMP_STACKSIZE:str = "OMP_STACKSIZE"
    """Environment variable corresponding to OMP stack size."""

    SUPPORT_EXTENSIONSS:List[str] = ".xyz .mol .sdf .pdb".split()
    """A List[str] of supported file extensions for xtb."""

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for xtb."""
        return [self.XTB_EXE]

    def required_env_vars(self)-> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return [self.KMP_STACKSIZE, self.OMP_STACKSIZE]


    def required_py_modules(self)-> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return list()


def default_xtb_config() -> XTBConfig:
    """Creates a deep-copied default version of the XTBConfig() class."""
    return deepcopy(XTBConfig())
