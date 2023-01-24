"""Specialization of the Residue() class for a Solvent. Has no additional attributes vs Residue() parent class.
Meant to be stored alongside other Residue() and Residue() derived objects (Ligand() and MetalUnit()) inside of 
the Chain() object. Solvent() objects SHOULD NOT exist on their own.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-05
"""
from __future__ import annotations
from copy import deepcopy

from .atom import Atom
from .residue import Residue
from typing import List
from enzy_htp.chemical import enum as renum


class Solvent(Residue):
    """Represents a specific Solvent in a .pdb file. Typically created from a base Residue() object using the
    residue_to_solvent() method found in enzy_htp.structure.solvent.py. Has no additional attributes vs
    the parent Residue() class. The value is_solvent() has been hardcoded to True and the the value
    for Solvent.rtype_ is set to ResidueType.SOLVENT. Meant to be stored alongside other Residue() and
    Residue()-derived classes (MetalUnit() and Ligand()) in Chain() objects.

    Attributes:
        (same as residue)
    """

    def __init__(self,
                 residue_idx: int,
                 residue_name: str,
                 atoms: List[Atom],
                 parent=None):
        """Constructor for Solvent. Identical to Residue() ctor. SHOULD NOT be called directly by users.
        Instead use enzy_htp.structure.solvent.residue_to_solvent()."""
        Residue.__init__(self, residue_idx, residue_name, atoms, parent)
        self.rtype = renum.ResidueType.SOLVENT

    def is_solvent(self) -> bool:
        """Checks if the Solvent() is an solvent. Hard-coded to True for this derived class."""
        return True

    # def clone(self) -> Solvent:
    #     """Creates deepcopy of self."""
    #     return deepcopy(self)

    def __str__(self) -> str:
        return f"Solvent({self._idx}, {self._name}, atom:{len(self._atoms)}, {self._parent})"


def residue_to_solvent(residue: Residue) -> Solvent:
    """Conversion method that creates a Solvent() instance from a ba6se Residue()."""
    return Solvent(residue.idx, residue.name, residue.atoms, residue.parent)
