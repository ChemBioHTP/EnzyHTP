"""Defines a BaseConfig() class that all config classes for supported software should inherit from. 
The Config() class found in enzy_htp/_config/config.py should NOT use this.
Class contains:
    
    + default constructor
    + required_executables() method
    + required_env_vars() method
    + required_py_modules() method
    + __getitem__ implementation
    + __setitem__ implementation

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-08-01
"""
import abc
from abc import ABC
from typing import Callable, Dict, List, Any

from enzy_htp import _LOGGER


class BaseConfig(ABC):
    """Abstract base calss for all Package Configs."""

    def __init__(self):
        """Trivial constructor"""
        pass
    
    @abc.abstractmethod
    def required_executables(self) -> List[str]:
        """A list of all required executables for <Package>Config."""
        pass

    @abc.abstractmethod
    def required_env_vars(self) -> List[str]:
        """A list of all required environment variables for <Package>Config."""
        pass

    @abc.abstractmethod
    def required_py_modules(self) -> List[str]:
        """A list of all required environment python modules for <Package>Config."""
        pass


    def __getitem__(self, key:str) -> Any:
        """Getter for the config object. Exits if the supplied key is invalid.

        Args:
            key: Name of the attribute as a str().

        Returns:
            The value associated with the attribute.
        """
        if not hasattr(self, key):
            _LOGGER.error(f"The Config object has not attribute {key}. Exiting...")
            exit( 1 )

        return getattr(self, key)

    def __setitem__(self, key:str, value: Any) -> None:
        """Setter for the config object. Checks if the attribute exists, exits if not.

        Args:
            key: Name of the attribute as a str().
            value: The value to set the attribute to.

        Returns:
            Nothing.
        """
        if not hasattr(self, key):
            _LOGGER.error(f"The Config object has not attribute {key}. Exiting...")
            exit( 1 )

        setattr(self, key, value )

