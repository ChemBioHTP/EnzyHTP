"""Definition for the Residue class. Residues are the most common unit of function within 
enzy_htp. A Residue() can be canonincal, non-canonical, solvent, or ligand. It is essentially
the catch all for PDB objects.
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
from __future__ import annotations
import copy
import sys
import math
from typing import Dict, Tuple, List, Union

import numpy as np
from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode
from enzy_htp.core import _LOGGER
from enzy_htp.core.exception import ResidueDontHaveAtom
import enzy_htp.chemical as chem
from enzy_htp.core.math_helper import get_geom_center

from .atom import Atom


class Residue(DoubleLinkedNode):
    """Most common functional unit in enzy_htp. Made up of Atom() objects and can be either
        canonical, non-canonical, solvent or a metal center.
    Attributes:
        (nessessary)
        children/atoms : A list of Atom() objects that make up the Residue().
        name : Residue name. (3-letter code)
        idx : The index of the Residue within the chain.
        parent/chain : Parent chain name.
        rtype : The ResidueType of the Residue().
    Derived properties:
        key(if_name) : the tuple key of residue's sequential identity e.g.: ("A", 1)
        num_atoms: number of belonging Atom()s
    """

    def __init__(self, residue_idx: int, residue_name: str, atoms: List[Atom], parent=None):
        """Constructor for the Residue() object"""
        self._name = residue_name
        self._idx = residue_idx
        self._rtype = chem.ResidueType.UNKNOWN
        self.set_parent(parent)
        self.set_children(atoms)

    #region === Getter-Attr ===
    @property
    def atoms(self) -> List[Atom]:
        """Returns a list of all Atom() objects that the Residue() "owns" """
        return self.get_children()

    @atoms.setter
    def atoms(self, val):
        self.set_children(val)

    @property
    def _atoms(self) -> List[Atom]:
        """alias for _children. Prevent change children but _atoms dont change"""
        return self._children

    @_atoms.setter
    def _atoms(self, val):
        self.set_children(val)

    @property
    def chain(self):
        """Getter for the Residue()'s parent chain id."""
        return self.get_parent()

    @chain.setter
    def chain(self, val):
        self.set_parent(val)

    @property
    def idx(self) -> int:
        """Getter for the Residue()'s index."""
        return self._idx

    @idx.setter
    def idx(self, val):
        self._idx = val

    @property
    def rtype(self) -> chem.ResidueType:
        """Getter for the Residue()'s chem.ResidueType value."""
        return self._rtype

    @rtype.setter
    def rtype(self, val):
        self._rtype = val

    @property
    def name(self) -> str:
        """Getter for the Residue()'s name. (3-letter code)"""
        return self._name

    @name.setter
    def name(self, val):
        self._name = val

    #endregion

    #region === Getter-Prop ===
    def key(self, if_name: bool = False) -> Tuple[str, int]:
        """
        an unique tuple to indicate a specific residue in the enzyme
        """
        if if_name:
            return (self.chain.name, self._idx, self._name)
        return (self.chain.name, self._idx)

    @property
    def key_str(self) -> str:
        """The Residue key as a str in the format <chain_name>.<residue_idx>. Note that if the given
        Residue does not have a parent chain, it is left blank (e.g. ".1"."""
        tokens:List[str] = ["", str(self.idx)]

        if self.chain is not None:
            tokens[0] = self.chain.name

        return ".".join(tokens)

    @property
    def seqres(self) -> chem.SeqRes:
        """Converts the Residue() into its respective SeqRes object."""
        return chem.SeqRes(
            model=1,
            chain=self.parent.name,
            idx=self.idx,
            seq_idx=None,
            name=self.name,
            missing=False
        )

    @property
    def num_atoms(self) -> int:
        """Number of atoms in the Residue."""
        return len(self._atoms)

    @property
    def sequence_name(self) -> str:
        """get the sequence name of the residue. 1-letter if canonical. 3-letter if non"""
        if self.is_canonical():
            return chem.convert_to_one_letter(self.name)
        return f" {self.name} "

    def has_atom_name(self, name:str) -> bool:
        """Does the Residue() contain an atom with the specified name? """
        for aa in self.atoms:
            if aa.name == name:
                return True

        return False
        

    def find_atom_name(self, name: str) -> Atom:
        """find child Atom base on its name"""
        result = list(filter(lambda a: a.name == name, self.atoms))
        if len(result) > 1:
            _LOGGER.error(f"residue {self} have more than 1 {name}")
            sys.exit(1)
        if len(result) == 0:
            raise ResidueDontHaveAtom(self, name, f"residue {self} dont have {name}")
        return result[0]

    @property
    def atom_name_list(self) -> List[str]:
        """get a list of atom names in the residue"""
        return list(map(lambda a: a.name, self.atoms))

    @property
    def atom_idx_list(self) -> List[int]:
        """get a list of atom indexes in the residue"""
        return list(map(lambda a: a.idx, self.atoms))

    @property
    def ca_coord(self) -> Tuple[float, float, float]:
        """get the C-Alpha coordinate of the residue"""
        atom_ca = self.find_atom_name("CA")
        return atom_ca.coord

    @property
    def geom_center(self) -> Tuple[float, float, float]:
        """get the geom_center coordinate of the residue"""
        return get_geom_center([i.coord for i in self.atoms])

    @property
    def element_composition(self) -> set:
        """get the element composition of the residue"""
        result = set()
        for atom in self.atoms:
            result.add(atom.element)
        return result

    @property
    def mainchain_atoms(self) -> List[Atom]:
        """return a list of mainchain atoms"""
        result = []
        atom_names = "C CA N".split() # TODO do we add CP1 CP2?
        for name in atom_names:
            result.append(self.find_atom_name(name))
        return result
    
    @property
    def atom_idx_mapper(self) -> Dict[int, Atom]:
        """the mapper for idx -> atom"""
        return {atom.idx : atom for atom in self.atoms}
    
    @property
    def hydrogens(self) -> List[Atom]:
        """Return all the hydrogen atoms in the Residue/Ligand."""
        result: List[Atom] = list(filter(lambda atom: atom.is_hydrogen() , self.atoms))
        return result

    def clash_count(self, other:Residue, radius:float=2.0, ignore_H:bool=True) -> int:
        """Counts how many clashes exist between the current residue and another. The distance cutoff 
        for each clash as well as whether Hydrogen atoms should be considered is determined by user parameters.

        Args:
            other: Ther other Residue() object to compare against.
            radius: How many Angstroms can separate two atoms before they are clashing? Default is 2.0
            ignore_H: Should Hydrogen atoms be considered for the clash count?

        Returns:
            The number of clashes between the Residue() and the other Residue().
        """
        if self == other:
            _LOGGER.warning("Attempted to count clashes between Residue and itself.")
            return 0

        count:int = 0
        for aa in other.atoms:
            if ignore_H and aa.element == 'H':
                continue
            else:
                count += self.clashes(aa, radius)

        return count
            
    def clashes(self, other:Atom, radius:float) -> bool:
        """Counts how many clashes exist between the current Residue and a given Atom(). The distance cutoff
        for each clash is determined by user supplied parameters.

        Args:
            other: The Atom() to check for clashes with.
            radius: How many Angstroms can separate two atoms before they are clashing?

        Returns:
            The number of clashes between the Atom() and the Residue().
        """
        for aa in self.atoms:
            if aa.distance_to(other) <= radius:
                return True
        return False

    def c_side_residue(self) -> Union["Residue", None]:
        """get the sibling Residue that connects to the C atom.
        return None if reached the chain terminal."""
        result = self.chain.find_residue_idx(self.idx + 1)
        if result is None:
            if self.chain.largest_res_idx > self.idx:
                _LOGGER.warning(f"{self} have no C-side residue but it is not a C-ter."
                                " Something wrong in your residue indexing. "
                                "It could be your index is not continous")
        return result

    def n_side_residue(self)-> Union["Residue", None]:
        """get the sibling Residue that connects to the N atom
        return None if reached the chain terminal."""
        result = self.chain.find_residue_idx(self.idx - 1)
        if result is None:
            if self.chain.smallest_res_idx < self.idx:
                _LOGGER.warning(f"{self} have no N-side residue but it is not a N-ter."
                                " Something wrong in your residue indexing. "
                                "It could be your index is not continous")
        return result
    
    def find_idxes_atom_list(self, atom_idx_list: int) -> List[Atom]:
        """find atom base on its idx. return a list reference of the atoms."""
        result = []
        atom_mapper = self.atom_idx_mapper
        for idx in atom_idx_list:
            result.append(atom_mapper[idx])
        return result
    #endregion

    #region === Checker ===
    def is_canonical(self) -> bool:
        """Checks if Residue() is canonical. Inherited by children."""
        return self._rtype == chem.ResidueType.CANONICAL

    def is_noncanonical(self) -> bool:
        """Checks if Residue() is canonical. Inherited by children."""
        return (self.is_modified() or self.is_ligand() or self.is_metal())

    def is_modified(self) -> bool:
        """Checks if Residue() is a modified residue. Inherited by children.
        TODO: would MODIFIED a better name?"""
        return self._rtype == chem.ResidueType.MODIFIED

    def is_ligand(self) -> bool:
        """Checks if Residue() is a ligand. Inherited by children."""
        return self._rtype == chem.ResidueType.LIGAND

    def is_solvent(self) -> bool:
        """Checks if the Residue() is an rd_sovlent as defined by enzy_htp.chemical.solvent.RD_SOLVENT"""
        return self._rtype == chem.ResidueType.SOLVENT

    def is_counterions(self, counterion_list: List[str]=None) -> bool:
        """Checks if the Residue() is 
        counter ion. default only Na+ and Cl- counters.
        This list can be overwritten by {counterion_list}"""
        if counterion_list is None:
            counterion_list = chem.COUNTER_ION_LIST
        return self.name in counterion_list

    def is_metal(self) -> bool:
        """Checks if Residue() is a metal. Inherited by children."""
        return self._rtype == chem.ResidueType.METAL

    def is_metal_center(self):
        """determine if current metal is a coordination center.
        TODO more consistant way to determine for 'boundary' metals like Mg2+"""
        return self.is_metal() and self.name in chem.METAL_CENTER_MAP # fix the bug that it think U as nucleic_acid is also metal center

    def is_trash(self) -> bool:
        """Checks if the Residue() is an rd_non_ligand as defined by enzy_htp.chemical.solvent.RD_NON_LIGAND_LIST"""
        return self._rtype == chem.ResidueType.TRASH

    def is_residue_cap(self) -> bool:
        """Checks if the Residue() is residue_cap"""
        return self._rtype == chem.ResidueType.RESIDUECAP

    def is_missing_chain(self) -> bool:
        """Whether the Residue() lacks a chain parent"""
        return self._parent is None

    def is_sequence_eq(self, other: Residue) -> bool:
        """Comparator that checks for sequence same-ness."""
        return self.key() == other.key()

    def is_deprotonatable(self) -> bool:
        """
        check if this residue can minus a proton in a pH range of 1-14.
        (ambiguous protonation state, potential deprotonation)
        Base on the residue.
        """
        return self.name in chem.residue.DEPROTONATION_MAPPER

    def is_hetatom_noproton(self) -> bool:
        """check if the residue contain no acidic proton on side chain hetero atoms"""
        return self.name in chem.residue.NOPROTON_LIST

    def is_connected(self) -> bool:
        """check if every atom in the residue have a connect record"""
        return math.prod([atom.is_connected() for atom in self.atoms])

    def has_init_charge(self) -> bool:
        """check if self has charge"""
        return math.prod([atom.has_init_charge() for atom in self.atoms])  
    
    def has_hydrogens(self) -> bool:
        """Does the residue contain hydrogen atoms?"""
        element_list = self.element_composition
        return ('H' in element_list)

    #endregion

    #region === Editor ===
    def sort_atoms(self):
        """
        sort children atoms with their atom idx
        sorted is always better than not but Residue() is being lazy here
        so only this function is called will it sorted
        """
        self._children.sort(key=lambda x: x.idx)

    def fix_atom_names(self):
        """
        Atom names should be unique in a residue to represent its connectivity.
        This method check the overall topology of the residue and assign atoms with
        valid names.
        For canonical AA, its should follow those CONNECTIVITY_MAPPER
        For non-canonical, every name should be at least unique
        """
        raise Exception  #TODO need to figure out how to determine the connectivity without the name

    def renumber_atoms(self, atom_idx_list: List[int]) -> None:  # TODO
        """Renumbers the Residue()'s Atom()'s to match atom_idx_list
        """
        idx = 0
        for atom in self._atoms:
            atom._idx = atom_idx_list[idx]
            idx += 1

    def remove_atoms_not_in_list(self, keep_name_list: List[str]):
        """remove atoms that do not in the name list of keeping.
        Make change in place"""
        ref_atom_list: List[Atom] = copy.copy(self.atoms)
        for atom in ref_atom_list:
            if atom.name not in keep_name_list:
                _LOGGER.debug(f"deleting {atom}")
                atom.delete_from_parent()

    def shift(self, point:List[float] ) -> None:
        """Shift all atoms by the supplied vector. Checks if the input is the correct
        size and throws if not.
        
        Args:
            point: An iterable, indexable object containing 3 floats().

        Returns:
            Nothing.
        """
        if len(point) != 3:
            raise TypeError(f"The supplied point is NOT of length 3 ({point})")

        for aa in self.atoms:
            orig = aa.coord
            updated = (orig[0]+point[0], orig[1]+point[1], orig[2]+point[2])
            aa.coord = updated

    #endregion


    #region === Special ===
    def __str__(self) -> str:
        return f"Residue({self._idx}, {self._name}, atom:{len(self._atoms)}, {self._parent})"

    def __repr__(self) -> str:  # TODO
        return str(self)

    #endregion
