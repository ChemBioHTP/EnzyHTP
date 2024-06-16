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
from .noncanonical_base import NonCanonicalBase
import enzy_htp.chemical as chem

from enzy_htp.core import _LOGGER


class Ligand(NonCanonicalBase):
    """Represents a specific Ligand found in a .pdb file. (#@shaoqz: a non-covalently binding small molecule to the protein part of the enzyme. decouple with PDB)
        Typically created from a base Residue() object using
        the residue_to_ligand() method found in enzy_htp.structure.ligand.py. In addition to base attributes, has
        net_charge attribute which is Union[float,None]. The value is_ligand() has been hard-coded to True and
        Ligand.rtype_ is set to ResidueType.LIGAND. Meant to be stored alongside other Residue() and Residue()-derived
        classes (MetalUnit() and Solvent()) in Chain() objects.

    Attributes:
        net_charge : The net charge of the molecule as an int.
        multiplicity: The multiplicity of the molecule as an int
        bonds: A List[Dict] containing bond information.
        conformer_coords: A List[List[Tuple[float,float,float]]] containing the coordinates of conformers.
    """

    def __init__(self,
                residue_idx: int,
                residue_name: str,
                atoms: List[Atom], 
                parent=None, 
                **kwargs):
        """
        Constructor for Ligand. Identical to Residue() ctor but also takes net_charge value.
        """
        NonCanonicalBase.__init__(self, residue_idx, residue_name, atoms, parent, **kwargs)
        self.rtype = chem.ResidueType.LIGAND

        self.bonds = kwargs.get('bonds', list())  # TODO change this to a more general represetation
        self.ligand_ensemble_ = kwargs.get('ligand_ensemble', None)
        # TODO: The idea for Structure and Ligand is to only have the math model recorded. 
        # But I see there is not really a better place to record the source of placement. 
        # Would a self.mimo attribute a good idea for these? So that we can distinguish those
        # that should not assumed to exist but just used in some protocols.
        # (similar to job_manager.ClusterJob)

    # === Getter-Attr ===

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
    
    @property
    def has_ensemble(self) -> bool:
        return self.ligand_ensemble_ is not None
    
    @property
    def ensemble(self) -> "LigandEnsemble":
        return self.ligand_ensemble_

    @ensemble.setter
    def ensemble(self, ensemble: "LigandEnsemble" ) -> None:
        self.ligand_ensemble_ = ensemble


    # === Special ===
    def __str__(self) -> str:
        return f"Ligand({self._idx}, {self._name}, atom:{len(self._atoms)}, {self._parent})"

def residue_to_ligand(residue: Residue, net_charge: float = None) -> Ligand:
    """Convenience function that converts Residue to ligand."""
    return Ligand(residue.idx, residue.name, residue.atoms, residue.parent, net_charge=net_charge)
