"""Defines a BaseInterface() class that all software interface classes should inherit from.
The Interface() class found in enzy_htp/_interface/interface.py should NOT use this. 
Class contains:
    
    + default constructor
    + getter for <Package>Config instance 
    + getter for whether the current environment is compatible with the <Package>Interface() 


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-06-07

"""

#TODO(CJ): need more documentation here

from abc import ABC
from typing import Callable, Dict, List

from enzy_htp.core import env_manager as em


class BaseInterface:
    """Abstract base class for all Package Interfaces. 

    Attributes:
        parent_ : The parent class of the <Package>Interface. Should ALWAYS be the main Interface.
        config_ : The relevant <Package>Config() class for the relevant Package.
        env_manager_ : The perfonal EnvironmentManager() class for the <Package>Interface.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config=None, default_config: Callable = None) -> None:
        """Simplistic constructor that optionally takes a <Package>Config object as its only argument.
        Also checks if the current environment is compatible with the <Package>Interface().
        """
        self.parent_ = parent
        self.config_ = config
        if not self.config_:
            self.config_ = default_config()
        self.env_manager_ = em.EnvironmentManager(env_vars=self.config_.required_env_vars(),
                                                  executables=self.config_.required_executables(),
                                                  py_modules=self.config_.required_py_modules())
        self.env_manager_.check_environment()
        self.compatible_env_ = self.env_manager_.is_missing()

    def parent(self):
        """Getter for the parent of the <Package>Interface."""
        return self.parent_

    def config(self):
        """Getter for the <Package>Config() instance belonging to the class."""
        return self.config_

    def compatible_environment(self) -> bool:
        """Checks if the current environment is compatible with all possible needs for the <Package>Interface.
        Returns:
                Whether the current environment is suitable for the <Package>Interface().
        """
        return self.compatible_env_

    def missing_executables(self) -> List[str]:
        """ """
        return self.env_manager_.missing_executables()

    def missing_env_vars(self) -> List[str]:
        """ """
        return self.env_manager_.missing_env_vars()

    def missing_py_modules(self) -> List[str]:
        """ """
        return self.env_manager_.missing_py_modules()
