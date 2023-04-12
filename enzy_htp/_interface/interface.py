"""Submodule describing an interface class which enables interfacing into all supported applications within EnzyHTP. This
model sees all interactions between EnzyHTP and a given application via a data attribute which is an instance of a 
specific application interface. At present, the below packages can be interacted with via the corresponding Interface
class. Packages:
    + AmberMD, AmberInterface
    + BCL, BCLInterface
    + Gaussian, GaussianInterface
    + MOE, MOEInterface
    + Multiwfn, MultiwfnInterface
    + PyMOL, PyMOLInterface
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-20
"""
from typing import List

from enzy_htp import core
from ..core.logger import _LOGGER
from .amber_interface import AmberInterface
from .bcl_interface import BCLInterface
from .gaussian_interface import GaussianInterface
from .moe_interface import MOEInterface
from .multiwfn_interface import MultiwfnInterface
from .pymol_interface import PyMOLInterface
from .rosetta_interface import RosettaInterface

from enzy_htp._config import Config


class Interface:
    """Wrapper class that houses access to individual <Package>Interface classes that exist in EnzyHTP.
    Each <Package>Interface is available as the attribute Interface.<package> (all lower case). Each instance
    needs an EnzyHTP.Config instance to hold all data.
    Attributes:
        amber: Corresponds to instance of AmberInterface().
        bcl: Corresponds to instance of BCLInterface().
        gaussian: Corresponds to instnce of GaussianInterface().
        moe: Corresponds to instance of MOEInterface().
        multiwfn: Corresponds to instance of MultiwfnInterface().
        pymol: Corresponds to instance of PyMOLInteface().
        rosetta: Corresponds to instance of RosettaInterface().
    """

    def __init__(self, config: Config):
        self.amber = AmberInterface(config._amber)
        self.bcl = BCLInterface(config._bcl)
        self.gaussian = GaussianInterface(config._gaussian)
        self.moe = MOEInterface(config._moe)
        self.multiwfn = MultiwfnInterface(config._multiwfn)
        self.pymol = PyMOLInterface(config._pymol)
        self.rosetta = RosettaInterface(config._rosetta)

        self.check_environment()


    def check_environment(self) -> None:
        """Checks for which elements are available in the environment. Gets executables and
        environment variables from children interface classes."""

        missing_exes:List[str] = list()
        missing_env_vars:List[str] = list()

        missing_exes.extend( self.amber.env_manager_.missing_executables() )
        missing_env_vars.extend( self.amber.env_manager_.missing_env_vars() )
        
        missing_exes.extend( self.bcl.env_manager_.missing_executables() )
        missing_env_vars.extend( self.bcl.env_manager_.missing_env_vars() )

        missing_exes.extend( self.gaussian.env_manager_.missing_executables() )
        missing_env_vars.extend( self.gaussian.env_manager_.missing_env_vars() )

        missing_exes.extend( self.moe.env_manager_.missing_executables() )
        missing_env_vars.extend( self.moe.env_manager_.missing_env_vars() )

        missing_exes.extend( self.multiwfn.env_manager_.missing_executables() )
        missing_env_vars.extend( self.multiwfn.env_manager_.missing_env_vars() )

        missing_exes.extend( self.pymol.env_manager_.missing_executables() )
        missing_env_vars.extend( self.pymol.env_manager_.missing_env_vars() )

        missing_exes.extend( self.rosetta.env_manager_.missing_executables() )
        missing_env_vars.extend( self.rosetta.env_manager_.missing_env_vars() )

        _LOGGER.info("Beginning environment check...")
        _LOGGER.info("Environment check complete!")
        
        missing:bool = False

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

        if missing:
            _LOGGER.warning("Some elements are missing in environment. Not all functionality will be possible.")
        else:
            _LOGGER.info("All elements are available in environment!")
