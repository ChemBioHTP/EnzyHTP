"""This module contains submodules that defines different "handle" types
that appear in EnzyHTP's interface of external software. The "handle" stands
for functions that solves a specific type of problem. For example, amber_interface
and openmm_interface will both have the MolDynStep handle defined here.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-19
"""

from .mol_dyn_step import MolDynStep
from .mol_dyn_parameterizer import MolDynParameterizer, MolDynParameter
from .mol_dyn_result import MolDynResult, MolDynResultEgg
from .qm_engine import QMOptimizationEngine, QMSinglePointEngine, QMResultEgg
from .modeling_engine import ModelingResultEgg, ModelingEngine
from .thermostability_engine import ddGFoldEngine, ddGResultEgg
