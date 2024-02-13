"""This module defines the interface of thermostability calculation engines as 
abstract classes.
Concrete classes of thermostability calculation engines may show up in interfacing 
modules of different external software.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-02-13
"""
from abc import ABC, abstractmethod
from typing import Union, List, Tuple

from ..base_interface import BaseInterface
from .modeling_engine import ModelingEngine, ModelingResultEgg

from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.structure import Structure
from enzy_htp.core.job_manager import ClusterJob

class ddGResultEgg(ModelingResultEgg):
    """This class defines the format of QM result eggs.
    These eggs are QM result place holders before the calculation.
    An example of eggs is file paths, so in Gaussian's case, it could be
    the .out file path etc. In general it should be some info to deduce
    where the result will be.
    The motivation of having this class is when running QM through ARMer
    the jobs.
    In most cases, the concrete class is a dataclass"""
    pass

class ddGFoldEngine(ModelingEngine):
    """the qm engine that performs single point calculation"""
    @abstractmethod
    def make_job(self, stru: Structure) -> Tuple[ClusterJob, ddGResultEgg]:
        """the method that make a ClusterJob that runs the QM"""
        pass

    @abstractmethod
    def run(self, stru: Structure) -> ElectronicStructure:
        """the method that runs the QM"""
        pass

    @abstractmethod
    def translate(self, egg: ddGResultEgg) -> ElectronicStructure:
        """the method convert engine specific results to general output"""
        pass
