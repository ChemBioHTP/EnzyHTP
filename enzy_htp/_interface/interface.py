"""Submodule describing an interface class which enables interfacing into all supported applications within EnzyHTP. This
model sees all interactions between EnzyHTP and a given application via a data attribute which is an instance of a 
specific application interface. At present, the below packages can be interacted with via the corresponding Interface
class. Packages:

    + AmberMD, AmberInterface
    + AlphaFill, AlphaFillInterface
    + BCL, BCLInterface
    + Gaussian, GaussianInterface
    + Modeller, ModellerInterface
    + MOE, MOEInterface
    + Mole2, Mole2Interface
    + Multiwfn, MultiwfnInterface
    + PyMol, PyMolInterface
    + RDKit, RDKitInterface
    + Rosetta, RosettaInterface
    + xtb, XTBInterface

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-20
"""
from typing import List, Dict

from enzy_htp import core
from ..core.logger import _LOGGER
from .amber_interface import amber_interface
from .alphafill_interface import AlphaFillInterface
from .bcl_interface import BCLInterface
from .gaussian_interface import gaussian_interface
from .modeller_interface import ModellerInterface
from .moe_interface import MOEInterface
from .mole2_interface import Mole2Interface
from .multiwfn_interface import MultiwfnInterface
from .pymol_interface import PyMolInterface
from .rdkit_interface import RDKitInterface
from .rosetta_interface import RosettaInterface
from .xtb_interface import xtb_interface 

from enzy_htp._config import Config


class Interface:
    """Wrapper class that houses access to individual <Package>Interface classes that exist in EnzyHTP.
    Each <Package>Interface is available as the attribute Interface.<package> (all lower case). Each instance
    needs an EnzyHTP.Config instance to hold all data.

    Attributes:
        _config: Corresponds to instance of Config().
        amber: Corresponds to instance of AmberInterface().
        alphafill: Corresponds to instance of AlphaFillInterface().
        bcl: Corresponds to instance of BCLInterface().
        gaussian: Corresponds to instance of GaussianInterface().
        modeller: Corresponds to instance of ModellerInterface().
        moe: Corresponds to instance of MOEInterface().
        mole2: Corresponds to instance of Mole2Interface().
        multiwfn: Corresponds to instance of MultiwfnInterface().
        pymol: Corresponds to instance of PyMolInteface().
        rdkit: Corresponds to instance of RDKitInterface().
        rosetta: Corresponds to instance of RosettaInterface().
        xtb: Corresponds to an isntance of XTBInterface().
    """

    def __init__(self, config: Config):
        """Constructor for the Interface(). Takes only a Config() class."""
        self._config = config
        self.alphafill = AlphaFillInterface(self, config._alphafill)
        amber_interface.set_parent(self) # made change for local access in amber_interface.py
        self.amber = amber_interface
        self.bcl = BCLInterface(self, config._bcl)
        gaussian_interface.set_parent(self)
        self.gaussian = gaussian_interface
        self.modeller = ModellerInterface(self, config._modeller)
        self.moe = MOEInterface(self, config._moe)
        self.mole2 = Mole2Interface(self, config._mole2)
        self.multiwfn = MultiwfnInterface(self, config._multiwfn)
        self.pymol = PyMolInterface(self, config._pymol)
        self.rdkit = RDKitInterface(self, config._rdkit)
        self.rosetta = RosettaInterface(self, config._rosetta)
        self.xtb = xtb_interface
        xtb_interface.set_parent(self)

        self.check_environment()
    
    @property
    def config(self) -> Config:
        """Getter for the config class."""
        return self._config

    def check_environment(self) -> None:
        """Checks for which elements are available in the environment. Gets executables and
        environment variables from children interface classes."""

        missing_exes: List[str] = list()
        missing_env_vars: List[str] = list()
        missing_py_modules: List[str] = list()

        for interface_name, ii in self.__dict__.items():
            if interface_name.startswith('_'):
                continue

            missing_exes.extend(ii.missing_executables())
            missing_env_vars.extend(ii.missing_env_vars())
            missing_py_modules.extend(ii.missing_py_modules())

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
