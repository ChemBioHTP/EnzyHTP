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

from .residue import Residue


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

    def insert_residue(self, new_res : Residue, sort_after : bool = True ) -> None:
        """Allows for insertion of a new Residue() object into the Chain. If the new Residue() is an exact
        copy, it fully overwrites the existing value. The sort_after flag specifies if the Residue()'s should
        be sorted by residue_number after the insertion.
        """
        for ridx, res in enumerate(self.residues_):
            if new_res.name == res.name and new_res.num_ == res.num_:
                self.residues_[idx] = deepcopy( new_res )
                break
        else:
            self.residues_.append( new_res )
        
        self.rename( self.name_ )
        
        if sort_after:
            self.residues_.sort( key=lambda r: r.num() )

    def is_metal(self) -> bool:
        """Checks if any metals are contained within the current chain."""
        return sum(list(map(lambda rr: rr.is_metal(), self.residues_)))

    def is_HET(self) -> bool:
        for rr in self.residues_:
            if not rr.is_canonical():
                return False
        return True

    def empty(self) -> bool:
        """Does the chain have any Residue()'s."""
        return len(self.residues_) == 0

    def residues(self) -> List[Residue]:
        """Access the child Residue() objects."""
        return self.residues_

    def name(self) -> str:
        """Getter for the Chain's name."""
        return self.name_

    def __getitem__(self, key: int) -> Residue:
        """Allows indexing into the child Residue() objects."""
        return self.residues_[key]

    def __delitem__(self, key: int) -> None:
        """Allows deleting of the child Residue() objects."""
        del self.residues_[key]

    def same_sequence(self, other: Chain) -> bool:
        """Comparison operator for use with other Chain() objects. Checks if residue list is identical in terms of residue name only."""
        self_residues: List[Residue] = self.residues_
        other_residues: List[Residue] = other.residues_
        # print(len(self_residues),'\t',len(other_residues))
        if len(self_residues) != len(other_residues):
            return False

        for s, o in zip(self_residues, other_residues):
            s: Residue
            o: Residue
            if not s.sequence_equivalent(o):
                return False
        return True

    def rename(self, new_name: str) -> None:
        """Renames the chain and propagates the new chain name to all child Residue()'s."""
        self.name_ = new_name
        res: Residue
        for ridx, res in enumerate(self.residues_):
            self.residues_[ridx].set_chain(new_name)

    def num_atoms(self) -> int:
        """Finds the total number of Atom() objects contained in the Residue() children objects."""
        total = 0
        for res in self.residues():
            total += res.num_atoms()
        return total

    def renumber_atoms(self, start: int = 1) -> int:
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
        for ridx, res in enumerate(self.residues_):
            idx = self.residues_[ridx].renumber_atoms(idx)
            idx += 1
        return idx - 1

    def get_pdb_lines(self) -> List[str]:
        """Generates a list of PDB lines for the Atom() objects inside the Chain(). Last line is a TER."""
        result = list()
        for res in self.residues_:
            result.extend(res.get_pdb_lines())
        result.append("TER")
        return result

    def __len__(self) -> int:
        """Returns number of Residue() or Residue()-dervied objects belonging to the Chain."""
        return len(self.residues_)


    def num_residues(self) -> int:
        """Returns number of Residue() or Residue()-dervied objects belonging to the Chain."""
        return len(self)

    def remove_residue(self, target_key : str ) -> None:
        """Given a target_key str of the Residue() residue_key ( "chain_id.residue_name.residue_number" ) format, 
        the Residue() is removed if it currently exists in the Chain() object."""
        for ridx, res in enumerate(self.residues_):
            if res.residue_key == target_key:
                break
        else:
            return
        
        del self.residues_[ridx]







