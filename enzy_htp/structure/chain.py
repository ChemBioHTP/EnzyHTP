"""Definition for the Chain class. Chains primarily store Residue() objects and organize them
within the overall structure of an enzyme.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-20
"""
from __future__ import annotations
from copy import deepcopy
from enzy_htp.core import _LOGGER
from typing import List

from enzy_htp.structure.atom import Atom

from .residue import Residue
#TODO(CJ): add a method for changing/accessing a specific residue

class Chain:
    """Class that represents a Chain of residues in a PDB file. Serves as a manager to the
    Residue() objects that it owns.

    Attributes:
        name_ : The name of the chain as a string.
        residues_ : A list of Residue() objects or derived types.
    """

    def __init__(self, name: str, residues: List[Residue]):
        """Initiation of a Chain with a name and list of residues."""
        self.name_ = name
        self.residues_: List[Residue] = deepcopy(residues)
        self.rename(self.name_)
   
    # === Getter-Attr (ref) ===
    def residues(self) -> List[Residue]:
        """Access the child Residue() objects."""
        return self.residues_

    def get_residue(self, traget_key: str) -> Residue:
        '''TODO: is there alt option for the key?'''
        pass

    def name(self) -> str:
        """Getter for the Chain's name."""
        return self.name_

    @property
    def atoms(self) -> List[Atom]:
        '''get all children Atoms'''
        result = list()
        for res in self:
            result.extend(res.atoms)
        return result

    # === Getter-Prop (cpy/new) ===
    def num_atoms(self) -> int:
        """Finds the total number of Atom() objects contained in the Residue() children objects."""
        total = 0
        for res in self.residues():
            total += res.num_atoms()
        return total

    def num_residues(self) -> int:
        """Returns number of Residue() or Residue()-dervied objects belonging to the Chain."""
        return len(self)

    # === Checker === 
    def is_metal(self) -> bool:
        """Checks if any metals are contained within the current chain."""
        return sum(list(map(lambda rr: rr.is_metal(), self.residues_)))

    def is_HET(self) -> bool:
        for rr in self.residues_:
            if not rr.is_canonical():
                return False
        return True #@shaoqz: why not use sum like above lol

    def is_empty(self) -> bool: #@shaoqz: @imp2 maybe name it is_empty?
        """Does the chain have any Residue()'s."""
        return len(self.residues_) == 0

    def is_same_sequence(self, other: Chain) -> bool:
        """Comparison operator for use with other Chain() objects. Checks if residue list is identical in terms of residue name only."""
        self_residues: List[Residue] = self.residues_
        other_residues: List[Residue] = other.residues_
        # print(len(self_residues),'\t',len(other_residues))
        if len(self_residues) != len(other_residues):
            return False

        for s, o in zip(self_residues, other_residues):
            s: Residue
            o: Residue
            if not s.is_sequence_equivalent(o): #@shaoqz: this is a good idea of having different levels of comparsion @imp2 after reading this method I found here it already comparing residues in the same position we only need to compare the name but not the key
                return False
        return True

    def is_same_coord(self, other: Chain) -> bool:
        '''check if self is same as other in coordinate of every atom'''
        self_atoms = self.atoms
        self_atoms.sort(key=lambda a: a.num)
        self_coord = map(lambda x:x.coord, self_atoms)
        other_atoms = other.atoms
        other_atoms.sort(key=lambda a: a.num)
        other_coord = map(lambda x:x.coord, other_atoms)
        for s, o in zip(self_coord, other_coord):
            if s != o:
                return False
        return True    
        
    # === Editor ===
    def add_residue(self, new_res: Residue, sort_after: bool = True, overwrite: bool = False) -> None: #@shaoqz: @imp overwriting a residue is rarely a demand but instead inserting one with same name and index but different stru is.
        """Allows for insertion of a new Residue() object into the Chain. If the new Residue() is an exact # TODO work on overwriting
        copy, it fully overwrites the existing value. The sort_after flag specifies if the Residue()'s should
        be sorted by residue_number after the insertion.
        """
        for ridx, res in enumerate(self.residues_):
            if new_res.name == res.name and new_res.num_ == res.num_: #@shaoqz: @imp should not work like this. Imaging the case where both residue is called LIG but they are different ligands. They are neither the same nor should be overwritten.
                self.residues_[ridx] = deepcopy(new_res)
                break
        else:
            self.residues_.append(new_res)

        self.rename(self.name_) #@shaoqz: maybe better to change this residue attribute only?

        if sort_after:
            self.residues_.sort(key=lambda r: r.num()) #@shaoqz: @imp this does not work when two different residue have the same index which often happens when you want to add a ligand to the structure.

    def remove_residue(self, target_key: str) -> None: #@shaoqz: @imp2 target key should not require the name. We should minimize the prerequisite of any input since its for HTP.
        """Given a target_key str of the Residue() residue_key ( "chain_id.residue_name.residue_number" ) format,
        the Residue() is removed if it currently exists in the Chain() object."""
        for ridx, res in enumerate(self.residues_):
            if res.residue_key == target_key:
                break
        else:
            return

        del self.residues_[ridx] #@shaoqz: why not move this line inside the loop?
    
    def rename(self, new_name: str) -> None: #@shaoqz: using the 2-way link sheet will get rid of functions like this but both works.
        """Renames the chain and propagates the new chain name to all child Residue()'s."""
        self.name_ = new_name
        res: Residue
        for ridx, res in enumerate(self.residues_):
            self.residues_[ridx].set_chain(new_name) #@shaoqz: why not just use res?

    def renumber_atoms(self, start: int = 1) -> int: #@shaoqz: @imp need to record the mapping of the index  #@shaoqz: also need one for residues
        """Renumbers the Atom()'s inside the chain beginning with "start" value and returns index of the last atom.
        Exits if start index <= 0.
        """
        if start <= 0:
            _LOGGER.error(
                f"Illegal start number '{start}'. Value must be >= 0. Exiting..."
            )
            exit(1)
        self.residues_ = sorted(self.residues_, key=lambda r: r.num())
        idx = start
        num_residues: int = self.num_residues()
        for ridx, res in enumerate(self.residues_):
            idx = self.residues_[ridx].renumber_atoms(idx)
            idx += 1
            terminal = (ridx < (num_residues - 1)) and (
                res.is_canonical() and not self.residues_[ridx + 1].is_canonical() #@shaoqz: @imp what does this mean? the TER line?
            )
            if terminal:
                idx += 1
        return idx - 1 

    # === Special ===
    def __getitem__(self, key: int) -> Residue:
        """Allows indexing into the child Residue() objects."""
        return self.residues_[key]

    def __delitem__(self, key: int) -> None:
        """Allows deleting of the child Residue() objects."""
        del self.residues_[key]

    def __len__(self) -> int:
        """Returns number of Residue() or Residue()-dervied objects belonging to the Chain."""
        return len(self.residues_)

    #region === TODO/TOMOVE ===
    def get_pdb_lines(self) -> List[str]: #@shaoqz: @imp move to the IO class
        """Generates a list of PDB lines for the Atom() objects inside the Chain(). Last line is a TER."""
        result = list()
        num_residues: int = self.num_residues()
        for idx, res in enumerate(self.residues_):
            terminal = (idx < (num_residues - 1)) and (
                res.is_canonical() and not self.residues_[idx + 1].is_canonical()
            )
            result.extend(res.get_pdb_lines(terminal)) #@shaoqz: why terminal?
        result.append("TER")
        return result
    #endregion
