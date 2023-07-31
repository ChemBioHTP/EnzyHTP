"""Submodule describing an interface class which enables interfacing into all supported applications within EnzyHTP. This
model sees all interactions between EnzyHTP and a given application via a data attribute which is an instance of a 
specific application interface. At present, the below packages can be interacted with via the corresponding Interface
class. Packages:
    + AmberMD, AmberInterface
    + BCL, BCLInterface
    + Gaussian, GaussianInterface
    + MOE, MOEInterface
    + Multiwfn, MultiwfnInterface
    + PyMol, PyMolInterface
    + Rosetta, RosettaInterface
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-20
"""
from typing import List, Dict

from enzy_htp import core
from ..core.logger import _LOGGER
from .amber_interface import AmberInterface
from .bcl_interface import BCLInterface
from .gaussian_interface import GaussianInterface
from .moe_interface import MOEInterface
from .multiwfn_interface import MultiwfnInterface
from .pymol_interface import PyMolInterface
from .rosetta_interface import RosettaInterface
from .pymol_interface import PyMolInterface

from enzy_htp._config import Config


class Interface:
    """Wrapper class that houses access to individual <Package>Interface classes that exist in EnzyHTP.
    Each <Package>Interface is available as the attribute Interface.<package> (all lower case). Each instance
    needs an EnzyHTP.Config instance to hold all data.
    Attributes:
        _config: TODO(CJ)
        amber: Corresponds to instance of AmberInterface().
        bcl: Corresponds to instance of BCLInterface().
        gaussian: Corresponds to instnce of GaussianInterface().
        moe: Corresponds to instance of MOEInterface().
        multiwfn: Corresponds to instance of MultiwfnInterface().
        pymol: Corresponds to instance of PyMolInteface().
        rosetta: Corresponds to instance of RosettaInterface().
    """

    def __init__(self, config: Config):
        self._config = config
        self.amber = AmberInterface(self, config._amber)
        self.bcl = BCLInterface(self, config._bcl)
        self.gaussian = GaussianInterface(self, config._gaussian)
        self.moe = MOEInterface(self, config._moe)
        self.multiwfn = MultiwfnInterface(self, config._multiwfn)
        self.pymol = PyMolInterface(self, config._pymol)
        self.rosetta = RosettaInterface(self, config._rosetta)

        self.check_environment()

    def config(self) -> Config:
        return self._config

    def check_environment(self) -> None:
        """Checks for which elements are available in the environment. Gets executables and
        environment variables from children interface classes."""

        #TODO(CJ): loop through __dict__ items

        missing_exes: List[str] = list()
        missing_env_vars: List[str] = list()
        missing_py_modules: List[str] = list()

        missing_exes.extend(self.amber.missing_executables())
        missing_env_vars.extend(self.amber.missing_env_vars())
        missing_py_modules.extend(self.amber.missing_py_modules())

        missing_exes.extend(self.bcl.missing_executables())
        missing_env_vars.extend(self.bcl.missing_env_vars())
        missing_py_modules.extend(self.bcl.missing_py_modules())

        missing_exes.extend(self.gaussian.missing_executables())
        missing_env_vars.extend(self.gaussian.missing_env_vars())
        missing_py_modules.extend(self.gaussian.missing_py_modules())

        missing_exes.extend(self.moe.missing_executables())
        missing_env_vars.extend(self.moe.missing_env_vars())
        missing_py_modules.extend(self.moe.missing_py_modules())

        missing_exes.extend(self.multiwfn.missing_executables())
        missing_env_vars.extend(self.multiwfn.missing_env_vars())
        missing_py_modules.extend(self.multiwfn.missing_py_modules())

        missing_exes.extend(self.pymol.missing_executables())
        missing_env_vars.extend(self.pymol.missing_env_vars())
        missing_py_modules.extend(self.pymol.missing_py_modules())

        missing_exes.extend(self.rosetta.missing_executables())
        missing_env_vars.extend(self.rosetta.missing_env_vars())
        missing_py_modules.extend(self.rosetta.missing_py_modules())

        _LOGGER.info("Beginning environment check...")
        _LOGGER.info("Environment check complete!")

        missing: bool = False

        if missing_exes:
            missing = True
            _LOGGER.warning("The following executables are missing:")
            for me in sorted(missing_exes):
                _LOGGER.warning(f"      {me}")

        if missing_env_vars:
            missing = True
            _LOGGER.warning("The following environment variables are missing:")
            for mev in sorted(missing_env_vars):
                _LOGGER.warning(f"      {mev}")

        if missing_py_modules:
            missing = True
            _LOGGER.warning("The following python modules are missing:")
            for mpm in sorted(missing_py_modules):
                _LOGGER.warning(f"      {mpm}")

        if missing:
            _LOGGER.warning("Some elements are missing in environment. Not all functionality will be possible.")
        else:
            _LOGGER.info("All elements are available in environment!")

    def check_for_elements(self, keys: List[str]) -> Dict[str, bool]:
        """ TODO"""
        #TODO(CJ)
        pass
