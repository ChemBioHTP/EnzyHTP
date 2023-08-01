"""Sub module that lays out the Config() class which stores all configuration settings
for the EnzyHTP package. Config() is meant to be used as a singleton and typically as the
enzy_htp.config attribute. DO NOT create instances of this object on your own. It is created
only by the module. For each <Package> that EnzyHTP interfaces with, there should be a
<Package>Interface in enzy_htp/_interface and a <Package>Config in enzy/_config. At present,
there are configuration settings for the below packages by the given <Package>Config:

    + AmberMD, AmberConfig
    + BCL, BCLConfig
    + Gaussian, GaussianConfig
    + MOE, MOEConfig
    + Multiwfn, MultiwfnConfig
    + PyMOL, PyMolConfig
    + Rosetta, RosettaConfig

In addition to specific packages, settings for the system are specified with SystemConfig

Note that for a given <Package>, there should exist a default_<package>_config function which
provides a default, deep-copied instance of the specific <Package>Confg.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-07-15
"""
from typing import Any, Dict
from .amber_config import AmberConfig, default_amber_config
from .bcl_config import BCLConfig, default_bcl_config
from .gaussian_config import GaussianConfig, default_gaussian_config
from .moe_config import MOEConfig, default_moe_config
from .multiwfn_config import MultiwfnConfig, default_multiwfn_config
from .pymol_config import PyMolConfig, default_pymol_config
from .rosetta_config import RosettaConfig, default_rosetta_config
from .system_config import SystemConfig, default_system_config

from enzy_htp.core import _LOGGER

# TODO(CJ): add more syntax checking for getter/setter
# TODO(CJ): need a BaseConfig() class so that the __getitem__ and __setitem__ can work better

class Config:
    """Class that holds all configuration settings for the different external softwares 
    that EnzyHTP relies on. Uses a parameter setting/getting syntax that uses [] operators
    as well as namespace dotting. Direct accession of the sub-configs is NOT recommended. 

    Attributes:
        _amber: Private instance of AmberConfig() with default settings.
        _bcl: Private instance of BCLConfig() with default settings.
        _gaussian: Private instance of GaussianConfig() with default settings.
        _moe: Private instance of MOEConfig() with default settings.
        _multiwfn: Private instance of MultiwfnConfig() with default settings.
        _pymol: Private instance of PyMolConfig() with default settings.
        _rosetta: Private instance of RosettaConfig() with default settings.
        _system: Private instance of SystemConfig() with default settings.
    """

    def __init__(self):
        """Constructor that creates a <Package>Config instance for each <package> using default_<package>_config."""
        self._amber = default_amber_config()
        self._bcl = default_bcl_config()
        self._gaussian = default_gaussian_config()
        self._multiwfn = default_multiwfn_config()
        self._moe = default_moe_config()
        self._pymol = default_pymol_config()
        self._rosetta = default_rosetta_config()
        self._system = default_system_config()

    def __getitem__(self, key: str) -> Any:
        """Getter for the settings in the Config() object. Uses the grammar: "<package>.<setting>" 
        Note that the <packge> is specified in all lowercase letters.

        Args:
            key: A string with the grammar "<package>.<name>" that you will returning.

        Returns:
            Corresponding value, if it exists.

        Raises:
            TypeError() if any part of the key is invalid.
        """
        #TODO(CJ): add the rosetta accessors
        if key.count("."):
            app, settings = key.split(".", 1)
            ptr = None
            if app == "amber":
                ptr = self._amber
            elif app == "bcl":
                ptr = self._bcl
            elif app == "gaussian":
                ptr = self._gaussian
            elif app == "multiwfn":
                ptr = self._multiwfn
            elif app == "pymol":
                ptr = self._pymol
            elif app == "rosetta":
                ptr = self._rosetta
            elif app == "system":
                ptr = self._system
            else:
                raise TypeError()
            return ptr[settings]
        else:
            raise TypeError()

    def __setitem__(self, key: str, value: Any) -> None:
        """Setter for the configuration values in the Config object. Uses the grammar
        config[key] = value, where key has the form "<package>.<name>" and value can be 
        any value. Note that the <packge> is specified in all lowercase letters.

        Args:
            key: A string with the grammar "<package>.<name>" for the value you will be updating.
            value: The value you are going to set the corresponding variable to. Can have any type.

        Raises:
            TypeError() if any part of the key is invalid.
        """
        #TODO(CJ): add the rosetta setters
        if key.count("."):
            app, settings = key.split(".", 1)
            ptr = None
            if app == "amber":
                ptr = self._amber
            elif app == "bcl":
                ptr = self._bcl
            elif app == "gaussian":
                ptr = self._gaussian
            elif app == "multiwfn":
                ptr = self._multiwfn
            elif app == "pymol":
                ptr = self._pymol
            elif app == "rosetta":
                ptr = self._rosetta
            elif app == "system":
                ptr = self._system
            else:
                raise TypeError()
            ptr[settings] = value
        else:
            raise TypeError()
