"""Specialization of the Residue() class for a Ligand. Primarly used to interface with PDB2PQR for ONLY the ligand
as it is removed from the full structure during protonation. In addition to Residue() object, has net_charge attribute.
Meant to be stored alongside other Residue() and Residue() derived objets (MetalUnit() and Solvent()) inside of the
Chain() object. Ligand() objects SHOULD NOT exist on their own. 

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
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


class Ligand(Residue):
    """Represents a specific Ligand found in a .pdb file. (#@shaoqz: a non-covalently binding small molecule to the protein part of the enzyme. decouple with PDB)
        Typically created from a base Residue() object using
        the residue_to_ligand() method found in enzy_htp.structure.ligand.py. In addition to base attributes, has
        net_charge attribute which is Union[float,None]. The value is_ligand() has been hard-coded to True and
        Ligand.rtype_ is set to ResidueType.LIGAND. Meant to be stored alongside other Residue() and Residue()-derived
        classes (MetalUnit() and Solvent()) in Chain() objects.

    Attributes:
        net_charge : The net charge of the molecule as an int.
        bonds: A List[Dict] containing bond information in the TRIPOS .mol2 format.
        conformer_coords: A List[List[Tuple[float,float,float]]] containing the coordinates of conformers.
    """

    def __init__(self, residue_idx: int, residue_name: str, atoms: List[Atom], parent=None, **kwargs):
        """
        Constructor for Ligand. Identical to Residue() ctor but also takes net_charge value.
        """
        self._net_charge = kwargs.get("net_charge", None)
        self._multiplicity = kwargs.get("multiplicity", None)
        self.bonds = kwargs.get('bonds', list()) 
        Residue.__init__(self, residue_idx, residue_name, atoms, parent)
        self.rtype = chem.ResidueType.LIGAND
        self.conformer_coords = list()

    
    def add_conformer(self, points:List[Tuple[float,float,float]] ) -> int:
        """Add the coordinates of a conformer to the Ligand(). Performs checks that the number of points matches the
        number of atoms in the Ligand. Also checks that each point is a tuple of size 3.
        
        Args:
            points: A List[Tuple[float,float,float]] containing the points of a new conformer.

        Returns:
            The number of conformers the Ligand has.
        """
        if len(points) != len(self.atoms):
            _LOGGER.error(f"A total of {len(points)} were supplied. Was expecting {len(self.atoms)}. Exiting...")
            exit( 1 )

        for row in points:
            if len(row) != 3:
                _LOGGER.error(f"The supplied point {row} is of the wrong dimension (expecting 3 elements). Exiting...")
                exit( 1 )

        self.conformer_coords.append( points )

        return self.n_conformers()


    def n_conformers(self) -> int:
        """How many conformers does this Ligand have?"""
        return len(self.conformer_coords) + 1

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

    # === Getter-Prop ===
    def clone(self) -> Ligand:
        """Creates deecopy of self."""
        return deepcopy(self)

    # === Checker ===
    def is_ligand(self) -> bool:
        """Checks if the Residue is a ligand. Always returns True for this specialization."""
        return True

    # === Editor ===
    def fix_atom_names(self) -> None:
        """
        Atom names should be unique in a ligand.
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
        return f"Ligand({self._idx}, {self._name}, atom:{len(self._atoms)}, {self._parent})"


def residue_to_ligand(residue: Residue, net_charge: float = None) -> Ligand:
    """Convenience function that converts Residue to ligand."""
    return Ligand(residue.idx, residue.name, residue.atoms, residue.parent, net_charge=net_charge)
