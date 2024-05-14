"""This module defines the interface of MolDynStep as an abstract class.
Concrete classes of MolDynStep may show up in interfacing module of different
MD software.
TODO could it be a concrete class too that handle some general operations?

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-19
"""

from abc import ABC, abstractmethod
from typing import List, Tuple, Union
from enzy_htp.core.job_manager import ClusterJob
from enzy_htp.structure.structure_constraint import StructureConstraint

from ..base_interface import BaseInterface
from .modeling_engine import ModelingEngine
from .mol_dyn_parameterizer import MolDynParameter
from .mol_dyn_result import MolDynResult, MolDynResultEgg

class MolDynStep(ModelingEngine):
    """A modular/indivisible step of Molecular Dynamics simulation.
    This is expected to be the most general building block of all different
    MD based sampling methods."""

    @property
    @abstractmethod
    def parent_interface(self) -> BaseInterface:
        """the interface of the corresponding engine software.
        normally contains all the constructors"""
        pass

    @property
    @abstractmethod
    def length(self) -> float:
        """the simulation length of this step (unit: ns)"""
        pass

    @property
    @abstractmethod
    def timestep(self) -> float:
        """the timestep of the simulation. (unit: ns)"""
        pass

    @property
    @abstractmethod
    def temperature(self) -> Union[float, List[Tuple[float]]]:
        """the temperature of the simulation. can be a list of 1d coordinates that indicate
        a changing temperature. e.g.: [(0,0), (0.5,300)] the 1st element is time in ns and
        the second is temperature at the time point."""
        pass

    @property
    @abstractmethod
    def thermostat(self) -> str:
        """the algorithm of the thermostat."""
        pass

    @property
    @abstractmethod
    def constrain(self) -> StructureConstraint:
        """the StructureConstraint object that indicates a geometry constrain in the step."""
        pass

    @property
    @abstractmethod
    def if_report(self) -> bool:
        """whether the step reports the output"""
        pass

    @abstractmethod
    def make_job(self, input_data: Union[MolDynParameter, MolDynResultEgg]) -> Tuple[ClusterJob, MolDynResultEgg]:
        """the method that make a ClusterJob that runs the step"""
        pass

    @abstractmethod
    def run(self, input_data: Union[MolDynParameter, MolDynResult]) -> MolDynResult:
        """the method that runs the step"""
        pass

    @abstractmethod
    def translate(self, egg: MolDynResultEgg) -> MolDynResult:
        """the method convert engine specific results to general output"""
        pass

    @classmethod
    @abstractmethod
    def try_merge_jobs(cls, job_list: List[ClusterJob]) -> List[ClusterJob]:
        """the classmethod that merge a list of jobs from MolDynStep to fewer jobs"""
        pass

