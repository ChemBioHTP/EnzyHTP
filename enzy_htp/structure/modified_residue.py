"""Specialization of the Residue() class for a Modified Residues. 
In addition to Residue() object, has net_charge etc. attributes.
Meant to be stored alongside other Residue() and Residue() derived objets (MetalUnit() and Solvent()) inside of the
Chain() object.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-12
"""
from __future__ import annotations

from copy import deepcopy

from enzy_htp.core.general import swapped_dict

from .atom import Atom
from typing import Dict, List
from .residue import Residue
from .noncanonical_base import NonCanonicalBase
import enzy_htp.chemical as chem
from enzy_htp.core import _LOGGER
import networkx as nx


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
        """return a list of mainchain atoms. If the mainchain is already set, return those atoms. Otherwise, 
           find the shortest path between N-term and C-term with find_mainchain(), then return it."""
        if len(self._mainchain_atoms) > 0:
            return self._mainchain_atoms
        else:
            self._mainchain_atoms = self.find_mainchain()
            return self._mainchain_atoms

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


    def find_mainchain(self) -> List[Atom]:
        """
        Finds the shortest path from the N-terminal to the C-terminal using connectivity information
        and breadth-first-search
        """

        graph = nx.Graph()
        graph.add_nodes_from([aa for aa in self.atoms])

        start_atom: Atom = self.find_atom_name("N")
        end_atom: Atom = self.find_atom_name("C")

        # create "graph" of atoms
        for aa in self.atoms:
            for oa in aa.connect:
                graph.add_edge(aa, oa[0])
        try:
            path = nx.shortest_path(graph, start_atom, end_atom)
        except nx.NodeNotFound:
            raise AttributeError(f"Path from n-term to c-term does not exist for {self.name}")
        
        return path
    
    def clone_connectivity(self, other):
        """clone connectivity from {other}.
        The atoms are aligned and a connect record is added to each atom referencing
        the connect from the aligned atom from other"""
        # 1. san check
        # - make sure other's atoms are not less than self's atoms

        if len(self.atom_name_list) > len(other.atom_name_list) or not set(self.atom_name_list).issubset(set(other.atom_name_list)):
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
                if other_cnt_atom in other_self_atom_mapper:
                    self_cnt_atom = other_self_atom_mapper[other_cnt_atom]
                    result.append((self_cnt_atom, bond_type))
            self_atom.connect = result
    
    # === Special ===
    def __str__(self) -> str:
        return f"ModifiedResidue({self._idx}, {self._name}, atom:{len(self._atoms)}, {self._parent})"


def residue_to_modified_residue(residue: Residue, net_charge: float = None) -> ModifiedResidue:
    """Convenience function that converts Residue to modified_residue."""
    return ModifiedResidue(residue.idx, residue.name, residue.atoms, residue.parent, net_charge=net_charge)
