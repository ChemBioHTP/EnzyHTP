"""This module defines the interface of MolDynParameterizer and MolDynParameter as abstract classes.
Concrete classes of MolDynParameterizer and MolDynParameter may show up in interfacing module of different
MD software.
TODO it could also contain some concrete methods

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-19
"""

from abc import ABC, abstractmethod

from enzy_htp.structure import Structure

class MolDynParameter(ABC):
    """The parameter of Molecular Dynamics simulation. Different package have
    parameters of different format which corresponds to different concrete classes"""

    @property
    @abstractmethod
    def engine(self) -> str:
        """the engine name that should be hardcoded in each concrete class"""
        pass


class MolDynParameterizer(ABC):
    """The parameterizer for Molecular Dynamics simulation.
    Solvation is also contained in this step"""

    @property
    @abstractmethod
    def engine(self) -> str:
        """the engine name that should be hardcoded in each concrete class"""
        pass

    @abstractmethod
    def run(self, stru: Structure) -> MolDynParameter:
        pass
