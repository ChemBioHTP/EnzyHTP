"""This module defines the interface of thermostability calculation engines as 
abstract classes.
Concrete classes of thermostability calculation engines may show up in interfacing 
modules of different external software.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-02-13
"""
from abc import abstractmethod
from typing import List, Tuple

from .modeling_engine import ModelingEngine, ModelingResultEgg

from enzy_htp.mutation_class import Mutation
from enzy_htp.structure import Structure
from enzy_htp.core.job_manager import ClusterJob

class ddGResultEgg(ModelingResultEgg):
    """see ModelingResultEgg"""
    pass

class ddGFoldEngine(ModelingEngine):
    """the ddg fold calculation engine that performs ddg fold calculation"""
    @abstractmethod
    def make_job(self, stru: Structure, mutant_space: List[Mutation]) -> Tuple[ClusterJob, ddGResultEgg]:
        """the method that make a ClusterJob that runs the ddg fold calculation"""
        pass

    @abstractmethod
    def run(self, stru: Structure) -> float:
        """the method that runs the ddg fold calculation"""
        pass

    @abstractmethod
    def translate(self, egg: ddGResultEgg) -> float:
        """the method convert engine specific results to general output"""
        pass
