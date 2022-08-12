"""Sub module that lays out the Config() class which stores all configuration settings
for the EnzyHTP package. Config() is meant to be used as a singleton and typically as the
enzy_htp.config attribute. DO NOT create instances of this object on your own. It is created
only by the module. For each <Package> that EnzyHTP interfaces with, there should be a
<Package>Interface in enzy_htp/_interface and a <Package>Config in enzy/_config. At present,
there are configuration settings for the below packages by the given <Package>Config:

    + AmberMD, AmberConfig
    + Gaussian, GaussianConfig
    + Multiwfn, MultiwfnConfig

In addition to specific packages, settings for the system are specified with SystemConfig

Note that for a given <Package>, there should exist a default_<package>_config function which
provides a default, deep-copied instance of the specific <Package>Confg.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-07-15
"""
from typing import Any, Dict

from .amber_config import AmberConfig, default_amber_config
from .gaussian_config import GaussianConfig, default_gaussian_config
from .multiwfn_config import MultiwfnConfig, default_multiwfn_config
from .system_config import SystemConfig, default_system_config

from enzy_htp.core import _LOGGER

# TODO(CJ): add more syntax checking for getter/setter


class Config:
    """Class that holds all configuration settings for the different external softwares 
    that EnzyHTP relies on. Uses a parameter setting/getting syntax that uses [] operators
    as well as namespace dotting. Direct accession of the sub-configs is NOT recommended. 

    Attributes:
        _amber: Private instance of AmberConfig() with default settings.
        _gaussian: Private instance of GaussianConfig() with default settings.
        _multiwfn: Private instance of MultiwfnConfig() with default settings.
        _system: Private instance of SystemConfig() with default settings.
    """

    def __init__(self):
        """Constructor that creates a <Package>Config instance for each <package> using default_<package>_config."""
        self._amber = default_amber_config()
        self._gaussian = default_gaussian_config()
        self._multiwfn = default_multiwfn_config()
        self._system = default_system_config()

    def __getitem__(self, key: str) -> Any:
        """Getter for the settings in the Config() object. Uses the grammar: "<package>.<setting>" 
        
        Args:
            key: A string with the grammar "<package>.<name>" that you will returning.

        Returns:
            Corresponding value, if it exists.

        Raises:
            TypeError() if any part of the key is invalid.
        """
        if key.count("."):
            app, settings = key.split(".", 1)
            ptr = None
            if app == "amber":
                ptr = self._amber
            elif app == "gaussian":
                ptr = self._gaussian
            elif app == "multiwfn":
                ptr = self._multiwfn
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
        any value:

        Args:
            key: A string with the grammar "<package>.<name>" for the value you will be updating.
            value: The value you are going to set the corresponding variable to. Can have any type.

        Raises:
            TypeError() if any part of the key is invalid.
        """
        if key.count("."):
            app, settings = key.split(".", 1)
            ptr = None
            if app == "amber":
                ptr = self._amber
            elif app == "gaussian":
                ptr = self._gaussian
            elif app == "multiwfn":
                ptr = self._multiwfn
            elif app == "system":
                ptr = self._system
            else:
                raise TypeError()
            ptr[settings] = value
        else:
            raise TypeError()
