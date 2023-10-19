"""Specialization of the Residue() class for a noncanonical parts (e.g.: ModifiedResidue or Ligand). 
In addition to Residue() object, has net_charge etc. attributes.
Meant to be stored alongside other Residue() and Residue() derived objets (MetalUnit() and Solvent()) inside of the
Chain() object.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-18
"""
from __future__ import annotations

from copy import deepcopy

import numpy as np

from .atom import Atom
from typing import List, Tuple
from .residue import Residue
import enzy_htp.chemical as chem

from enzy_htp.core import file_system as fs
from enzy_htp.core import (
    UnsupportedFileType,
    _LOGGER,
)


class NonCanonicalBase(Residue):
    """Base class for residue units other that canonical amino acid and water.
    Children classes include Ligand, ModifiedResidue, MetalUnit.
    Should never be instantiated.

    Additional Attributes:
        net_charge : The net charge of the molecule as an int.
        multiplicity: The multiplicity of the molecule as an int
    """

    def __init__(self, residue_idx: int, residue_name: str, atoms: List[Atom], parent=None, **kwargs):
        """
        Constructor for Ligand. Identical to Residue() ctor but also takes net_charge value.
        """
        self._net_charge = kwargs.get("net_charge", None)
        self._multiplicity = kwargs.get("multiplicity", None)
        Residue.__init__(self, residue_idx, residue_name, atoms, parent)
        self.rtype = chem.ResidueType.LIGAND

    # === Getter-Attr ===
    @property
    def net_charge(self) -> int:
        """Getter for the net_charge attribute."""
        return self._net_charge

    @net_charge.setter
    def net_charge(self, val: int):
        """Setter for the net_charge attribute."""
        self._net_charge = val

    @property
    def multiplicity(self) -> int:
        """Getter for the multiplicity attribute."""
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, val: int):
        """Setter for the multiplicity attribute."""
        self._multiplicity = val

    # === Editor ===
    def clone_connectivity(self, other: NonCanonicalBase) -> None:
        """clone connectivity from {other}."""
        pass #TODO
