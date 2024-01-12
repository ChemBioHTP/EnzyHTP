"""This module defines classes that describe the level of theory.
The level of theory means the theory model used in the modeling task.
Example: b3lyp-d3/def2-svp/SMD(water) is a level of theory using DFT method.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-28"""
from dataclasses import dataclass
from typing import Union
from abc import ABC, abstractmethod

@dataclass
class LevelOfTheory(ABC):
    """the base class of level of theory.
    These classes define the level of theory under
    a given theory framework such as QM or MM."""
    @abstractmethod
    def lot_type(self) -> str:
        return None

@dataclass
class QMLevelOfTheory(LevelOfTheory):
    """the level of theory of a QM simulation.
    basis_set: 
        the basis function set.
    method: 
        the method for solving the Schrödinger equation. (e.g.: HF, wb97xd, ...)
    solvent: 
        the name of the solvent or a dictionary that specific the parameters
        of the solvent. It needs to be consistent with the solvation model.
    solv_method:
        the name of the solvation model."""
    basis_set: str
    method: str
    solvent: Union[str, dict] = None
    solv_method: str = None

    def lot_type(self) -> str:
        return "qm"

@dataclass
class MMLevelOfTheory(LevelOfTheory):
    """the level of theory of a MM based simulation.
    force_field: 
        the name of the force field
    ligand_method: 
        the method parameterizating the ligand"""
    force_field: str
    ligand_method: str

    def lot_type(self) -> str:
        return "mm"
