"""Specialization of the Residue() class for a Solvent. Has no additional attributes vs Residue() parent class.
Meant to be stored alongside other Residue() and Residue() derived objects (Ligand() and MetalAtom()) insdie of 
the Chain() object. Solvent() objects SHOULD NOT exist on their own.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
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
    the parent Residue() class. The value is_rd_solvent() has been hardcoded to True and the the value
    for Solvent.rtype_ is set to ResidueType.SOLVENT. Meant to be stored alongside other Residue() and
    Residue()-derived classes (MetalAtom() and Ligand()) in Chain() objects.

    Attributes:

    """

    def __init__(self, residue_key: str, atoms: List[Atom]):
        """Constructor for Solvent. Identical to Residue() ctor. SHOULD NOT be called directly by users.
        Instead use enzy_htp.structure.solvent.residue_to_solvent()."""
        Residue.__init__(self, residue_key, atoms)
        self.set_rtype(renum.ResidueType.SOLVENT)

    def is_rd_solvent(self) -> bool:
        """Checks if the Solvent() is an rd solvent. Hard-coded to True for this derived class."""
        return True

    def clone(self) -> Solvent:
        """Creates deepcopy of self."""
        return deepcopy(self)


def residue_to_solvent(ptr: Residue) -> Solvent:
    """Conversion method that creates a deepcopied Solvent() instance from a base Residue()."""
    return deepcopy(Solvent(ptr.residue_key, ptr._atoms))
