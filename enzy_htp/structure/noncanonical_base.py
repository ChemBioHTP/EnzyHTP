"""Specialization of the Residue() class for a noncanonical parts (e.g.: ModifiedResidue or Ligand). 
In addition to Residue() object, has net_charge etc. attributes.
Meant to be stored alongside other Residue() and Residue() derived objets (MetalUnit() and Solvent()) inside of the
Chain() object.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-18
"""
from __future__ import annotations

from typing import List
from .residue import Residue
from .atom import Atom

import enzy_htp.chemical as chem
from enzy_htp.core.general import swapped_dict
from enzy_htp.core import _LOGGER


class NonCanonicalBase(Residue):
    """Base class for Residue units other that canonical amino acid and water.
    Parent class: Residue
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
        """clone connectivity from {other}.
        The atoms are aligned and a connect record is added to each atom referencing
        the connect from the aligned atom from other"""
        # 1. san check
        # - make sure self and other have same set of atom names
        if set(self.atom_name_list) != set(other.atom_name_list):
            _LOGGER.error(f"Atom names are not consistent between {self} and {other}. Cannot align.")
            raise ValueError
        # - make sure other is connected
        if not other.is_connected():
            _LOGGER.error(f"clone target not connected: {other}.")
            raise ValueError

        # 2. align
        self_other_atom_mapper = {}
        for atom in self.atoms:
            self_other_atom_mapper[atom] = other.find_atom_name(atom.name)
        other_self_atom_mapper = swapped_dict(self_other_atom_mapper)

        # 3. clone
        for self_atom in self.atoms:
            result = []
            other_atom: Atom = self_other_atom_mapper[self_atom]
            for other_cnt_atom, bond_type in other_atom.connect:
                self_cnt_atom = other_self_atom_mapper[other_cnt_atom]
                result.append((self_cnt_atom, bond_type))
            self_atom.connect = result
            