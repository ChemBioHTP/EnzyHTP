"""Definition for the Chain class. Chains primarily store Residue() objects and organize them
within the overall structure of an enzyme.
Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-20
"""
from __future__ import annotations
from copy import deepcopy
import itertools
import sys
from enzy_htp.core import _LOGGER
from typing import Iterable, List, Tuple, Union
from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode

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
            _LOGGER.warning(f"residue idx {idx} out of chain's range {self}"
                            )  #TODO may be make this an error
            return None
        return result[0]

    #endregion

    #region === Getter-Prop ===
    def residue_idx_interval(self,
                             if_str: bool = True
                             ) -> Union[str, Iterable[Tuple[int, int]]]:
        """
        a range representation of containing residue indexes
        Args & Returns:
            if_str:
                0: return a iterable of Tuples [(1,20),(25,30)]
                1: return a str "1-20,25-30"
        """
        interval_list = get_interval_from_list(self.residue_idxs)
        if if_str:
            range_strs = map(lambda x: f"{x[0]}-{x[1]}", interval_list)
            return ",".join(range_strs)
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
    def chain_type(self) -> str:
        """
        Returns:
            return the chain type. Possible chain types are:
                peptide
                Any composition from metal, ligand, solvent, trash
        """
        chain_type = []
        if self.is_peptide():
            return "peptide"
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
    def is_peptide(self) -> bool:
        """
        if there is any residue not canonical
        """
        return not sum(
            list(
                map(lambda rr: (not rr.is_canonical()) and
                    (not rr.is_noncanonical()), self._residues)))

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
            if not s.is_sequence_eq(
                    o
            ):  #@shaoqz: this is a good idea of having different levels of comparsion @imp2 after reading this method I found here it already comparing residues in the same position we only need to compare the name but not the key
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

    #endregion

    # === Editor ===
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

    def add_residue(
        self,
        new_res: Residue,
        sort_after: bool = True,
        overwrite: bool = False
    ) -> None:  # TODO#@shaoqz: @imp overwriting a residue is rarely a demand but instead inserting one with same name and index but different stru is.
        """Allows for insertion of a new Residue() object into the Chain. If the new Residue() is an exact # TODO work on overwriting
        copy, it fully overwrites the existing value. The sort_after flag specifies if the Residue()"s should
        be sorted by residue_number after the insertion.
        """
        for ridx, res in enumerate(self._residues):
            if new_res.name == res.name and new_res.num_ == res.num_:  #@shaoqz: @imp should not work like this. Imaging the case where both residue is called LIG but they are different ligands. They are neither the same nor should be overwritten.
                self._residues[ridx] = deepcopy(new_res)
                break
        else:
            self._residues.append(new_res)

        self.rename(
            self._name)  #@shaoqz: maybe better to change this residue attribute only?

        if sort_after:
            self._residues.sort(
                key=lambda r: r.idx()
            )  #@shaoqz: @imp this does not work when two different residue have the same index which often happens when you want to add a ligand to the structure.

    def remove_residue(
        self, target_key: str
    ) -> None:  # TODO#@shaoqz: @imp2 target key should not require the name. We should minimize the prerequisite of any input since its for HTP.
        """Given a target_key str of the Residue() residue_key ( "chain_id.residue_name.residue_number" ) format,
        the Residue() is removed if it currently exists in the Chain() object."""
        for ridx, res in enumerate(self._residues):
            if res.residue_key == target_key:
                break
        else:
            return

        del self._residues[ridx]  #@shaoqz: why not move this line inside the loop?

    def rename(
        self, new_name: str
    ) -> None:  # TODO#@shaoqz: using the 2-way link sheet will get rid of functions like this but both works.
        """Renames the chain and propagates the new chain name to all child Residue()"s."""
        self._name = new_name
        res: Residue
        for ridx, res in enumerate(self._residues):
            self._residues[ridx].set_chain(new_name)  #@shaoqz: why not just use res?

    def renumber_atoms(
        self,
        start: int = 1
    ) -> int:  # TODO#@shaoqz: @imp need to record the mapping of the index  #@shaoqz: also need one for residues
        """Renumbers the Atom()"s inside the chain beginning with "start" value and returns index of the last atom.
        Exits if start index <= 0.
        """
        if start <= 0:
            _LOGGER.error(
                f"Illegal start number '{start}'. Value must be >= 0. Exiting...")
            exit(1)
        self._residues = sorted(self._residues, key=lambda r: r.idx())
        idx = start
        num_residues: int = self.num_residues
        for ridx, res in enumerate(self._residues):
            idx = self._residues[ridx].renumber_atoms(idx)
            idx += 1
            terminal = (ridx < (num_residues - 1)) and (
                res.is_canonical() and not self._residues[ridx + 1].is_canonical(
                )  #@shaoqz: @imp what does this mean? the TER line?
            )
            if terminal:
                idx += 1
        return idx - 1

    #region === Special ===
    def __str__(self):
        """
        concise string representation of the chain
        """
        return f"Chain({self._name}, residue: {self.residue_idx_interval()})"

    #endregion


# TODO go to core
def get_interval_from_list(target_list: List[int]) -> Iterable[Tuple[int, int]]:
    """
    convert a list of int to the interval/range representation
    Returns:
        a generater of tuples with each indicating the start/end of the interval
    Example:
        >>> list(get_interval_from_list([1,2,3,6,7,8]))
        [(1,3),(6,8)]
    reference: https://stackoverflow.com/questions/4628333
    """
    # clean input
    target_list = sorted(set(target_list))
    # here use enum id as a ref sequence and group by the deviation
    for i, j in itertools.groupby(
            enumerate(target_list),
            lambda ref_vs_target: ref_vs_target[1] - ref_vs_target[0]):
        j = list(j)
        yield j[0][1], j[-1][1]
