"""Definition for the Chain class. Chains primarily store Residue() objects and organize them
within the overall structure of an enzyme.
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com> 
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-20
"""
from __future__ import annotations
from copy import deepcopy
import math
import sys
from typing import Iterable, List, Tuple, Union

from enzy_htp.core import _LOGGER
from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode
from enzy_htp.core.general import get_interval_from_list

from enzy_htp.structure.atom import Atom
import enzy_htp.chemical as chem

from .residue import Residue


class Chain(DoubleLinkedNode):
    """Class that represents a Chain of residues in a PDB file.
    Attributes:
        name : The name of the chain as a string.
        children/residues : A list of Residue() objects or derived types.
        parent/protein : the parent protein
    Derived properties:
        atoms : composing atoms
        num_atoms
        num_residues
        chain_type
        residue_idx_interval : the containing residue index range
    """

    def __init__(self, name: str, residues: List[Residue], parent=None):
        """Initiation of a Chain with a name and list of residues."""
        self._name = name
        self.set_parent(parent)
        self.set_children(residues)

    #region === Getter-Attr ===
    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, val: str):
        self._name = val

    @property
    def protein(self):
        return self.get_parent()

    @protein.setter
    def protein(self, val):
        self.set_parent(val)

    @property
    def residues(self) -> List[Residue]:
        """Access the child Residue() objects."""
        return self.get_children()

    @residues.setter
    def residues(self, val: List[Residue]):
        self.set_children(val)

    @property
    def _residues(self) -> List[Residue]:
        """alias for _children. prevent changing _children but _residues holds the same"""
        return self._children

    @_residues.setter
    def _residues(self, val: List[Residue]):
        self.set_children(val)

    def find_residue_idx(self, idx: int) -> Union[Residue, None]:
        """find the residue with correponding index"""
        result = list(filter(lambda r: r.idx == idx, self.children))
        if len(result) > 1:
            _LOGGER.error(f"chain {self} have more than 1 residue on {idx}")
            sys.exit(1)
        if len(result) == 0:
            _LOGGER.warning(f"residue idx {idx} out of chain's range {self}")  #TODO may be make this an error
            return None
        return result[0]

    #endregion

    #region === Getter-Prop ===
    def residue_idx_interval(self, if_str: bool = True) -> Union[str, Iterable[Tuple[int, int]]]:
        """
        a range representation of containing residue indexes
        Args & Returns:
            if_str:
                0: return a iterable of Tuples [(1,20),(25,30)]
                1: return a str "1-20,25-30"
        """
        interval_list = get_interval_from_list(self.residue_idxs)
        if if_str:
            contain_list = [f"{x[0]}-{x[1]}" if x[0] != x[1] else f"{x[0]}" for x in interval_list]
            if len(contain_list) == 1:
                return contain_list[0]
            range_strs = ",".join(contain_list)
            return range_strs
        return interval_list

    @property
    def atoms(self) -> List[Atom]:
        """get all children Atoms"""
        result = []
        for res in self:
            result.extend(res.atoms)
        return result

    @property
    def num_atoms(self) -> int:
        """Finds the total number of Atom() objects contained in the Residue() children objects."""
        total = 0
        for res in self._residues:
            total += res.num_atoms
        return total

    @property
    def num_residues(self) -> int:
        """Returns number of Residue() or Residue()-dervied objects belonging to the Chain."""
        return len(self)

    @property
    def residue_idxs(self) -> List[int]:
        """return a list of indexes of containing residues"""
        return list(map(lambda x: x.idx, self._residues))

    @property
    def largest_res_idx(self) -> int:
        """return the largest residue index in the chain"""
        return max(self.residue_idxs)

    @property
    def smallest_res_idx(self) -> int:
        """return the smallest residue index in the chain"""
        return min(self.residue_idxs)

    def c_ter_residue(self) -> Residue:
        """get the C-ter residue of the self chain.
        TODO this needs to be more general same as nter"""
        return self.residues[-1]

    def n_ter_residue(self) -> Residue:
        """get the N-ter residue of the self chain"""
        return self.residues[0]

    @property
    def chain_type(self) -> str:
        """
        Returns:
            return the chain type. Possible chain types are:
                peptide
                Any composition from metal, ligand, solvent, trash
        """
        chain_type = []
        if self.is_polypeptide():
            return "polypeptide"
        if self.has_metal():
            chain_type.append("metal")
        if self.has_ligand():
            chain_type.append("ligand")
        if self.has_trash():
            chain_type.append("trash")
        if self.has_solvent():
            chain_type.append("solvent")
        return ",".join(chain_type)

    @property
    def sequence(self) -> str:
        """return the sequence of current chain in string"""
        result = ""
        self.sort_residues()
        res: Residue
        for res in self._residues:
            result += res.sequence_name
        return result.strip()

    #endregion

    #region === Checker ===
    def is_polypeptide(self) -> bool:
        """
        if there is any non-aminoacid part in chain
        """
        return not sum(list(map(lambda rr: (not rr.is_canonical()) and (not rr.is_modified()), self._residues)))

    def has_metal(self) -> bool:
        """Checks if any metals are contained within the current chain."""
        return sum(list(map(lambda rr: rr.is_metal(), self._residues)))

    def has_ligand(self) -> bool:
        return sum(list(map(lambda rr: rr.is_ligand(), self._residues)))

    def has_solvent(self) -> bool:
        return sum(list(map(lambda rr: rr.is_solvent(), self._residues)))

    def has_trash(self) -> bool:
        return sum(list(map(lambda rr: rr.is_trash(), self._residues)))

    def is_HET(self) -> bool:
        for rr in self._residues:
            if not rr.is_canonical():
                return False
        return True

    def is_empty(self) -> bool:
        """Does the chain have any Residue()"s."""
        return len(self._residues) == 0

    def is_same_sequence(self, other: Chain) -> bool:  # TODO
        """Comparison operator for use with other Chain() objects. Checks if residue list is identical in terms of residue name only."""
        self_residues: List[Residue] = self._residues
        other_residues: List[Residue] = other.residues
        # print(len(self_residues),"\t",len(other_residues))
        if len(self_residues) != len(other_residues):
            return False

        for s, o in zip(self_residues, other_residues):
            s: Residue
            o: Residue
            if not s.name == o.name:
                return False
        return True

    def is_same_coord(self, other: Chain, tol: float = 0.01) -> bool:
        """Checks if this Chain() and another Chain() have the same atoms and respective coordinates.
        Args:
            other: A different Chain() object to compare to.
            tol: Numeric tolerance for difference in each xyz value. Default is 0.01.
        Returns:
            If the Atom()'s owned by each Chain() are equivalent.
        """
        self_atoms = self.atoms
        self_atoms.sort(key=lambda a: a.idx)
        self_coord = map(lambda x: x.coord, self_atoms)
        other_atoms = other.atoms
        other_atoms.sort(key=lambda a: a.idx)
        other_coord = map(lambda x: x.coord, other_atoms)
        for s, o in zip(self_coord, other_coord):
            if abs(s[0] - o[0]) > tol or abs(s[1] - o[1]) > tol or abs(s[2] - o[2]) > tol:
                return False
        return True

    def is_connected(self) -> bool:
        """check whether all atoms within the chain have connected initiated"""
        return math.prod([atom.is_connected() for atom in self.atoms])

    def has_init_charge(self) -> bool:
        """check if self has charge"""
        return math.prod([atom.has_init_charge() for atom in self.atoms])  
    #endregion

    #region === Editor ===
    def remove_trash(self):
        """
        remove trash ligands in the chain
        """
        for i in range(len(self._residues) - 1, -1, -1):
            if self._residues[i].rtype == chem.ResidueType.TRASH:
                _LOGGER.debug(f"removing TRASH: {self._residues[i]}")  # pylint: disable=logging-fstring-interpolation
                del self._residues[i]

    def sort_residues(self):
        """
        sort children residues with their residue idx
        sorted is always better than not but Chain() is being lazy here
        so only when sort is nessessary will it sort
        """
        self._children.sort(key=lambda x: x.idx)

    def add(self, new_res: Residue,
            overwrite: bool = False, sort: bool = True,) -> None:
        """Method that adds a new residue based on the idx. If idx is None it is added in the end by default.
        Args:
            overwrite: if overwrite when residue same idx. if not overwrite the new
                res will be inserted after the idx (which will cause the idx change of all following residues).
            sort: if sort after adding"""
        new_res.parent = self
        push_number = 0
        if (new_res.idx is None) or (new_res.idx not in self.residue_idxs):
            self.residues.append(new_res)
        else:
            for i, res in enumerate(self.residues):
                if res.idx == new_res.idx:
                    if overwrite:
                        self.residues[i] = new_res
                    else:
                        self.residues.insert(i, new_res)
                        push_number = 1
                    break

        if push_number:
            push_idx = new_res.idx
            for res in self.residues:
                if res.idx > push_idx:
                    res.idx += 1
            new_res.idx += 1

        if sort:
            self.sort_residues()
    #endregion

    #region === Special ===
    def __str__(self):
        """
        concise string representation of the chain
        """
        return f"Chain({self._name}, residue: {self.residue_idx_interval()})"

    #endregion
