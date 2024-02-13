"""This module defines the base interface of modeling engine as 
an abstract class.
Concrete classes of modeling engine are for example MolDynStep, QMSinglePointEngine.
etc.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-02-13
"""
from abc import ABC, abstractmethod
from typing import Tuple, Any

from enzy_htp.structure import Structure
from enzy_htp.core.job_manager import ClusterJob

class ModelingResultEgg(ABC):
    """This class defines the format of modeling result eggs.
    These eggs are modeling result place holders before the calculation.
    An example of eggs is file paths, so in Gaussian's case, it could be
    the .out file path etc. In general it should be some info to deduce
    where the result will be.
    The motivation of having this class is when running modeling through ARMer
    and the submission and execution of the jobs are handled outside of the 
    engine.
    In most cases, the concrete class is a dataclass"""
    pass

class ModelingEngine(ABC):
    """the engine of modeling in general"""
    @property
    @abstractmethod
    def engine(self) -> str:
        """the engine name that should be hardcoded in each concrete class"""
        pass

    @abstractmethod
    def make_job(self, stru: Structure) -> Tuple[ClusterJob, ModelingResultEgg]:
        """the method that make a ClusterJob that runs the modeling"""
        pass

    @abstractmethod
    def run(self, stru: Structure) -> Any:
        """the method that runs the modeling locally.
        Should return a general result for example ElectronicStructure
        in QM's case."""
        pass

    @abstractmethod
    def translate(self, egg: ModelingResultEgg) -> Any:
        """the method convert engine specific result eggs to general result.
        For example, QMSinglePointEgg -> ElectronicStructure"""
        pass

