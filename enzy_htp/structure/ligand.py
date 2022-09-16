"""Specialization of the Residue() class for a Ligand. Primarly used to interface with PDB2PQR for ONLY the ligand
as it is removed from the full structure during protonation. In addition to Residue() object, has net_charge attribute.
Meant to be stored alongside other Residue() and Residue() derived objets (MetalUnit() and Solvent()) inside of the
Chain() object. Ligand() objects SHOULD NOT exist on their own.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""
from __future__ import annotations


from copy import deepcopy

import numpy as np

from .atom import Atom
from typing import List
from .residue import Residue
import enzy_htp.chemical as chem

from enzy_htp.core import file_system as fs
from enzy_htp.core import (
    UnsupportedFileType,
    _LOGGER,
)


class Ligand(Residue):
    """Represents a specific Ligand found in a .pdb file. (#@shaoqz: a non-covalently binding small molecule to the protein part of the enzyme. decouple with PDB)
        Typically created from a base Residue() object using
        the residue_to_ligand() method found in enzy_htp.structure.ligand.py. In addition to base attributes, has
        net_charge attribute which is Union[float,None]. The value is_ligand() has been hard-coded to True and
        Ligand.rtype_ is set to ResidueType.LIGAND. Meant to be stored alongside other Residue() and Residue()-derived
        classes (MetalUnit() and Solvent()) in Chain() objects.

    Attributes:
        net_charge : The net charge of the molecule as an int.
    """

    def __init__(self, residue_idx: int, residue_name: str, atoms: List[Atom], parent=None, **kwargs):
        """
        Constructor for Ligand. Identical to Residue() ctor but also takes net_charge value.
        """
        self.net_charge = kwargs.get("net_charge", None)
        Residue.__init__(self, residue_idx, residue_name, atoms, parent)
        self.rtype = chem.ResidueType.LIGAND

    # === Getter-Attr (ref) ===
    # === Getter-Prop (cpy/new) ===
    def get_net_charge(self) -> int:
        """Getter for the net_charge attribute."""
        return self.net_charge

    def clone(self) -> Ligand:
        """Creates deecopy of self."""
        return deepcopy(self)

    # === Checker === 
    def is_ligand(self) -> bool:
        """Checks if the Residue is a ligand. Always returns True for this specialization."""
        return True

    # === Editor ===
    def set_residue_number(self, num: int) -> None:
        """Changes the resdiue number for all of the substituent atoms."""  #@shaoqz: also set for itself?
        for idx, aa in enumerate(self._atoms):
            self._atoms[idx].residue_number = num

    # === Special ===

    #region === TODO/TOMOVE ===
    def build(self, out_path: str) -> None: #@shaoqz: to IO ; also it should be the same as for residue
        """Method that builds the given ligand to the specified path, making sure it is a pdb filepath."""
        ext = fs.get_file_ext(out_path).lower()
        if ext != ".pdb":
            _LOGGER.error(
                f"The supplied file path '{out_path}' does not ahve a '.pdb' extension. Exiting..."
            )
            exit(1)
        lines = list(
            map(
                lambda pr: pr[1].to_pdb_line(a_id=pr[0] + 1, c_id=" "),
                enumerate(self._atoms),
            )
        ) + ["TER", "END"]
        fs.write_lines(out_path, lines)
    #endregion


def residue_to_ligand(residue: Residue, net_charge: float = None) -> Ligand:
    """Convenience function that converts Residue to ligand."""
    return Ligand(residue.idx, residue.name, residue.atoms, residue.parent, net_charge=net_charge)
