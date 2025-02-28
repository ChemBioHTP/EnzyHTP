"""Specialization of the Residue() class for a Modified Residues. 
In addition to Residue() object, has net_charge etc. attributes.
Meant to be stored alongside other Residue() and Residue() derived objets (MetalUnit() and Solvent()) inside of the
Chain() object.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-12
"""
from __future__ import annotations

from copy import deepcopy

from .atom import Atom
from typing import Dict, List
from .residue import Residue
from .noncanonical_base import NonCanonicalBase
import enzy_htp.chemical as chem
from enzy_htp.core import _LOGGER


class ModifiedResidue(NonCanonicalBase):
    """Represents a specific ModifiedResidue found in a .pdb file. (#@shaoqz: a non-covalently binding small molecule to the protein part of the enzyme. decouple with PDB)
        Typically created from a base Residue() object using
        the residue_to_modified_residue() method found in enzy_htp.structure.modified_residue.py. In addition to base attributes, has
        net_charge attribute which is Union[float,None]. The value is_modified_residue() has been hard-coded to True and
        ModifiedResidue.rtype_ is set to ResidueType.MODIFIED. Meant to be stored alongside other Residue() and Residue()-derived
        classes (MetalUnit() and Solvent()) in Chain() objects.

    Attributes:
        net_charge : The net charge of the molecule as an int.
    """

    def __init__(self, residue_idx: int, residue_name: str, atoms: List[Atom], parent=None, **kwargs):
        """
        Constructor for ModifiedResidue. Identical to Residue() ctor but also takes net_charge value.
        """
        NonCanonicalBase.__init__(self, residue_idx, residue_name, atoms, parent, **kwargs)
        self.rtype = chem.ResidueType.MODIFIED

        # mainchian
        mainchain_atom_names = kwargs.get("mainchain_atom_names", list())
        self.set_mainchain_atoms(mainchain_atom_names)

    # === Getter-Attr ===
    # === Getter-Prop ===
    def clone(self) -> ModifiedResidue:
        """Creates deecopy of self."""
        return deepcopy(self)

    @property
    def mainchain_atoms(self) -> List[Atom]:
        """return a list of mainchain atoms"""
        if len(self._mainchain_atoms) > 0:
            return self._mainchain_atoms
        else:
            raise AttributeError(f"Missing mainchain information for {self.name} at {self.key_str}! "
                                  "You can assign it using assign_mod_aa_mainchain()")

    # === Checker ===
    def is_modified_residue(self) -> bool:
        """Checks if the Residue is a modified_residue. Always returns True for this specialization."""
        return True

    # === Editor ===
    def set_mainchain_atoms(self, mainchain_atom_names: List[str]) -> List[Atom]:
        """set the mainchain_atoms by a list of names"""
        result = []
        for name in mainchain_atom_names:
            result.append(self.find_atom_name(name))
        self._mainchain_atoms = result

    def fix_atom_names(self) -> None:
        """
        Atom names should be unique in a modified_residue.
        This method assign atoms with valid names.
        """
        name_list = self.atom_name_list
        new_name_list = chem.get_valid_generic_atom_name(name_list)
        for name, atom in zip(new_name_list, self.atoms):
            if atom.name != name:
                _LOGGER.info(f"found atom with invalid name {atom}. changing it to {name}")
                atom.name = name

    # === Special ===
    def __str__(self) -> str:
        return f"ModifiedResidue({self._idx}, {self._name}, atom:{len(self._atoms)}, {self._parent})"


def residue_to_modified_residue(residue: Residue, net_charge: float = None) -> ModifiedResidue:
    """Convenience function that converts Residue to modified_residue."""
    return ModifiedResidue(residue.idx, residue.name, residue.atoms, residue.parent, net_charge=net_charge)
