"""the module for interface definition of strutcure file I/O classes

Author: shaoqz, <shaoqz@icloud.com>
Date: 2022-08-30
"""

from abc import ABC, abstractmethod
from typing import Union
from ..structure import Structure


class StructureParserInterface(ABC):
    """
    classes inherit from this class should have a unified interface for
    parsing file in/out Structure objects *note this is for file only*
    """

    @abstractmethod
    def __init__(self) -> None:
        """
        prevent it from instantiation
        """
        pass

    @classmethod
    @abstractmethod
    def get_structure(cls, path: str) -> Union[Structure, tuple]:
        """
        interface for converting file into the Structure().

        (* If it should be classmethod is doubltful but i use
        it here to keep consistance with the job manager)
        """
        pass

    @classmethod
    @abstractmethod
    def get_file_str(cls, stru: Structure) -> str:
        """
        interface for converting Structure() into the file content.

        (* If it should be classmethod is doubltful but i use
         it here to keep consistance with the job manager
         * file writing is not expected to be placed here as it should
         be handle by other modules)
        """
        pass
