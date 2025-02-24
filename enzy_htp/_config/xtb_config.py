"""Defines a  XTBConfig() which holds configuration settings for enzy_htp to interface with the xtb
software package. File also contains default_bcl_config() which creates a default version of 
the XTBConfig() class.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-08-17
"""
#TODO(CJ): update this stuff
from typing import List, Union
from copy import deepcopy

from .base_config import BaseConfig


class XTBConfig(BaseConfig):
    """Class that holds default alues for running xtb within enzy_htp.

    Attributes:
        XTB_EXE : str() corresponding to xtb application.
        N_ITER : int() corresponding to number of scf iterations to run.
    """

    XTB_EXE: str = "xtb"
    """Name of the xtb executable as a str()."""

    N_ITER: int = 500
    """Number of SCF iterations to try as an int()."""

    N_PROC: int = 1
    """Number of parrallel processes to use during a calculation as an int()."""

    KMP_STACKSIZE: str = "KMP_STACKSIZE"
    """Environment variable corresponding to KMP stack size."""

    OMP_STACKSIZE: str = "OMP_STACKSIZE"
    """Environment variable corresponding to OMP stack size."""

    SUPPORTED_EXTENSIONS: List[str] = ".xyz .mol .sdf .pdb".split()
    """A List[str] of supported file extensions for xtb."""
    
    SOLVATION_METHODS:List[str] = "ALPB GBSA".split()
    """A List[str] of supported solvation methods."""

    ALPB_SOLVENTS:List[str] = "acetone acetonitrile aniline benzaldehyde benzene ch2cl2 chcl3 cs2 dioxane dmf dmso ether ethylacetate furane hexandecane hexane methanol nitromethane octanol woctanol phenol toluene thf water".split()
    """A List[str] of supported solvents for the ALPB solvation method."""

    GBSA_SOLVENTS:List[str] = " acetone, acetonitrile, benzene, CH2Cl2, CHCl3, CS2, DMF (only GFN2-xTB), DMSO, ether, H2O, methanol, n-hexane THF and toluene".split()
    """A List[str] of supported solvents for the GBSA solvation method. Note that not all solvents will work with various xtb theory levels."""

    SUPPORTED_XTB_THEORY_LEVELS:List[str] = "GFN0 GFN1 GFN2".split()
    """Allowed levels of theory for xtb to use."""

    FORCE_CONSTANT:Union[None,float]=None
    #TODO(CJ)

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for xtb."""
        return [self.XTB_EXE]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return [self.KMP_STACKSIZE, self.OMP_STACKSIZE]

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required enviornment variables for xtb."""
        return list()


def default_xtb_config() -> XTBConfig:
    """Creates a deep-copied default version of the XTBConfig() class."""
    return deepcopy(XTBConfig())
