"""Definition for the Residue class. Residues are the most common unit of function within 
enzy_htp. A Residue() can be canonincal, non-canonical, solvent, or ligand. It is essentially
the catch all for PDB objects.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
from __future__ import annotations
import numpy as np
from plum import dispatch
from copy import deepcopy
from typing import Tuple, List

from .atom import Atom
from enzy_htp.core import _LOGGER
import enzy_htp.chemical as chem


class Residue:
    """Most common functional unit in enzy_htp. Made up of Atom() objects and can be either
	canonical, non-canonical, solvent or a metal center. 

    Attributes:
        self.atoms : A list of Atom() objects that make up the Residue().
        self.residue_key : Unique string identifier with format "chain.name.num".
		self.chain_ : Parent chain name.
		self.name : Residue name.
		self.num_ : The index of the Residue within the chain.
		self.rtype_ : The ResidueType of the Residue().
		self.min_line_ : The lowest one-indexed line of the children atoms.
		self.max_line_ : The highest one-indexed line of the children atoms.
    """

    def __init__(self, residue_key: str, atoms: List[Atom]):
        """Constructor for the Residue() object. Takes residue_key in format of "chain_name.resdue_name.residue_id" and a list of Atom() objects."""
        self.atoms = atoms
        self.residue_key = residue_key
        (chain, name, num) = self.residue_key.split(".")
        self.chain_ = chain
        self.name = name
        self.num_ = int(num)
        self.rtype_ = chem.ResidueType.UNKNOWN
        line_idxs = np.array(list(map(lambda aa: aa.line_idx, self.atoms)))

        if line_idxs.any():
            self.min_line_ = np.min(line_idxs)
            self.max_line_ = np.max(line_idxs)

    def renumber_atoms(self, start: int = 1) -> int:
        """Renumbers the Residue()'s Atom()'s beginning with "start" paramter, defaulted to 1. Returns the index of the last Atom().
		NOTE: errors if "start" is <= 0.
		"""
        if start <= 0:
            _LOGGER.error(
                f"Illegal start number '{start}'. Value must be >= 0. Exiting..."
            )
            exit(1)
        aa: Atom
        self.atoms = sorted(self.atoms, key=lambda aa: aa.atom_number)

        for idx, aa in enumerate(self.atoms):
            self.atoms[idx].atom_number = idx + start
        return idx + start

    def num_atoms(self) -> int:
        """Number of atoms in the Residue."""
        return len(self.atoms)

    def atom_list(self) -> List[Atom]:
        """Returns a list of all Atom() objects that the Residue() "owns" """
        return self.atoms

    def empty_chain(self) -> bool:
        """Whether the Residue() lacks a chain identitifer. NOT a chain object accessor."""
        return not len(self.chain_.strip())

    def chain(self) -> str:
        """Getter for the Residue()'s parent chain id."""
        return self.chain_

    def set_chain(self, val: str) -> None:
        """Sets chain and all child atoms to new value. Re-generates the self.residue_key attribute."""
        chain = val
        for idx in range(len(self.atoms)):
            self.atoms[idx].chain_id = val
        self.chain_ = val
        self.residue_key = f"{self.chain_}.{self.name}.{self.num_}"

    def num(self) -> int:
        """Getter for the Residue()'s number."""
        return self.num_

    def min_line(self) -> int:
        """The lowest one-indexed line one of the Residue()'s atoms come from."""
        return self.min_line_

    def max_line(self) -> int:
        """The highest one-indexed line one of the Residue()'s atoms come from."""
        return self.max_line_

    def line_range(self) -> Tuple[int, int]:
        """A tuple of the he lowest and highest one-indexed lines that the Residue()'s atoms come from."""
        return (self.min_line(), self.max_line())

    def neighbors(self, other: Residue) -> bool:
        """Checks if two residues from the same PDB are neighbors."""
        return (abs(self.min_line() - other.max_line()) == 1) or (
            abs(self.max_line() - other.min_line()) == 1
        )

    def is_ligand(self) -> bool:
        """Checks if Residue() is a ligand. Inherited by children."""
        return False

    def is_metal(self) -> bool:
        """Checks if Residue() is a metal. Inherited by children."""
        return False

    def is_canonical(self) -> bool:
        """Checks if Residue() is canonical. Inherited by children."""
        return self.name in chem.THREE_LETTER_AA_MAPPER

    def is_rd_solvent(self) -> bool:
        """Checks if the Residue() is an rd_sovlent as defined by enzy_htp.chemical.solvent.RD_SOLVENT"""
        return self.name in chem.RD_SOLVENT_LIST

    def is_rd_non_ligand(self) -> bool:
        """Checks if the Residue() is an rd_non_ligand as defined by enzy_htp.chemical.solvent.RD_NON_LIGAND_LIST"""
        return self.name in chem.RD_NON_LIGAND_LIST

    def sequence_equivalent(self, other: Residue) -> bool:
        """Comparator that checks for sequence same-ness."""
        return self.residue_key == other.residue_key

    def set_rtype(self, new_rtype: chem.ResidueType) -> None:
        """Setter for the Residue()'s chem.ResidueType value."""
        self.rtype_ = new_rtype

    def rtype(self) -> ResidueType:
        """Getter for the Residue()'s chem.ResidueType value."""
        return self.rtype_

    def get_pdb_lines(self) -> List[str]:
        """Method that gets the PDB lines for all of the Atom() objects in the Residue() object."""
        result = list()
        for aa in self.atoms:
            result.append(aa.to_pdb_line())
        return result

    def __str__(self) -> str:
        """String representationt that just shows residue key in format "chain_id.residue_name.residue_num" """
        return self.residue_key

    def __repr__(self) -> str:
        """String representationt that just shows residue key in format "chain_id.residue_name.residue_num" """
        return str(self)

    def clone(self) -> Residue:
        """Creates a deepcopy of self."""
        return deepcopy( self )
