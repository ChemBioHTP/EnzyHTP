"""

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-08-01
"""

from abc import ABC
from typing import Callable, Dict, List, Any


class BaseConfig(ABC):
    """ """

    
    @abstractmethod
    def required_executables(self) -> List[str]:
        """ """
        pass

    @abstractmethod
    def required_env_vars(self) -> List[str]:
        """ """
        pass

    @abstractmethod
    def required_py_modules(self) -> List[str]:
        """ """
        pass


    def __getitem__(self, key:str) -> Any:
        """ """
        return getattr(self, key)

    def __setitem__(self, key:str, value: Any) -> None:
        """ """
        setattr(self, key, value )

