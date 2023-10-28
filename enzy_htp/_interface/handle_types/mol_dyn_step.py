"""This module defines the interface of MolDynStep as an abstract class.
Concrete classes of MolDynStep may show up in interfacing module of different
MD software.
TODO could it be a concrete class too that handle some general operations?

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-19
"""

from abc import ABC, abstractmethod
from typing import List
from enzy_htp.core.job_manager import ClusterJob
from enzy_htp.core.mol_dyn_result import MolDynResult

from ..base_interface import BaseInterface

class MolDynStep(ABC):
    """A modular/indivisible step of Molecular Dynamics simulation.
    This is expected to be the most general building block of all different
    MD based sampling methods."""

    @property
    @abstractmethod
    def engine(self) -> str:
        """the engine name that should be hardcoded in each concrete class"""
        pass

    @property
    @abstractmethod
    def parent_interface(self) -> BaseInterface:
        """the interface of the corresponding engine software.
        normally contains all the constructors"""
        pass

    @property
    @abstractmethod
    def if_report(self) -> bool:
        """whether the step reports the output"""
        pass

    @abstractmethod
    def make_job(self) -> ClusterJob:
        """the method that make a ClusterJob that runs the step"""
        pass

    @abstractmethod
    def translate(self) -> MolDynResult:
        """the method convert engine specific results to general output"""
        pass

    @classmethod
    @abstractmethod
    def try_merge_jobs(cls, job_list: List[ClusterJob]) -> List[ClusterJob]:
        """the classmethod that merge a list of jobs from MolDynStep to fewer jobs"""
        pass

