"""This module defines the interface of MolDynParameterizer and MolDynParameter as abstract classes.
Concrete classes of MolDynParameterizer and MolDynParameter may show up in interfacing module of different
MD software.
TODO it could also contain some concrete methods

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-19
"""
from typing import Callable, List
from abc import ABC, abstractmethod

from enzy_htp.structure import Structure
from ..base_interface import BaseInterface

class MolDynParameter(ABC):
    """The parameter of Molecular Dynamics simulation. Different package have
    parameters of different format which corresponds to different concrete classes"""

    @property
    @abstractmethod
    def engine(self) -> str:
        """the engine name that should be hardcoded in each concrete class"""
        pass

    @property
    @abstractmethod
    def topology_file(self) -> str:
        """return the path of the topology file that composes the parameter"""
        pass

    @property
    def topology_parser(self) -> Callable:
        """return the parser object for topology_file"""
        pass

    @property
    def file_list(self) -> List[str]:
        """return a list of files that composes the parameter"""
        pass

class MolDynParameterizer(ABC):
    """The parameterizer for Molecular Dynamics simulation.
    Solvation is also contained in this step"""

    @property
    @abstractmethod
    def engine(self) -> str:
        """the engine name that should be hardcoded in each concrete class"""
        pass

    @property
    @abstractmethod
    def parameterizer_temp_dir(self) -> str:
        """the directory that contain the generated parm files"""
        pass

    @parameterizer_temp_dir.setter
    @abstractmethod
    def parameterizer_temp_dir(self, value: str) -> None:
        """set the directory that contain the generated parm files"""
        pass

    @property
    @abstractmethod
    def parent_interface(self) -> BaseInterface:
        """the interface of the corresponding engine software.
        normally contains all the constructors"""
        pass

    @abstractmethod
    def run(self, stru: Structure) -> MolDynParameter:
        """the parameterization process"""
        pass
