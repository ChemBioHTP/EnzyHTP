"""This module defines classes that describe the level of theory.
The level of theory means the theory model used in the modeling task.
Example: b3lyp-d3/def2-svp/SMD(water) is a level of theory using DFT method.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-28"""
from dataclasses import dataclass

@dataclass
class LevelofTheory:
    def lot_type(self) -> str:
        return None

@dataclass
class QMLevelofTheory(LevelofTheory):
    basis_set: str
    method: str
    solvent: str
    solv_method: str

    def lot_type(self) -> str:
        return "qm"

@dataclass
class MMLevelofTheory(LevelofTheory):
    force_field: str
    ligand_method: str

    def lot_type(self) -> str:
        return "mm"
