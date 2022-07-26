"""Definition for the Residue class. Residues are the most common unit of function within 
enzy_htp. A Residue() can be canonincal, non-canonical, solvent, or ligand. It is essentially
the catch all for PDB objects.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
# TODO(CJ): figure out how to inherit docstrings to the children classes.
from __future__ import annotations
import numpy as np
from collections import defaultdict
from plum import dispatch
from copy import deepcopy
from typing import Tuple, List, Dict

from .atom import Atom
from enzy_htp.core import _LOGGER
import enzy_htp.chemical as chem


class Residue:
    """Most common functional unit in enzy_htp. Made up of Atom() objects and can be either
        canonical, non-canonical, solvent or a metal center.

    Attributes:
        atoms_ : A list of Atom() objects that make up the Residue().
        residue_key : Unique string identifier with format "chain.name.num".
        chain_ : Parent chain name.
        name : Residue name.
        num_ : The index of the Residue within the chain.
        rtype_ : The ResidueType of the Residue().
        min_line_ : The lowest one-indexed line of the children atoms.
        max_line_ : The highest one-indexed line of the children atoms.
    """

    def __init__(self, residue_key: str, atoms: List[Atom]):
        """Constructor for the Residue() object. Takes residue_key in format of "chain_name.resdue_name.residue_id" and a list of Atom() objects."""
        self.atoms_ = atoms
        self.residue_key = residue_key #@shaoqz: @imp to identify a residue only chain_id and residue_id inside of the chain is enough. The residue name seem redundant especially when using it to aquire a residue. After reading more code i see you are using the key for 1) select the residue 2) comparing the residue 3) identifying the residue
        (chain, name, num) = self.residue_key.split(".")
        self.chain_ = chain
        self.name = name
        self.num_ = int(num)
        self.rtype_ = chem.ResidueType.UNKNOWN
        line_idxs = np.array(list(map(lambda aa: aa.line_idx, self.atoms_)))

        if line_idxs.any():
            self.min_line_ = np.min(line_idxs)
            self.max_line_ = np.max(line_idxs)

    # === Getter-Attr (ref) ===
    @property
    def atoms(self) -> List[Atom]:
        """Returns a list of all Atom() objects that the Residue() "owns" """
        return self.atoms_

    def chain(self) -> str:
        """Getter for the Residue()'s parent chain id."""
        return self.chain_

    def num(self) -> int:
        """Getter for the Residue()'s number."""
        return self.num_

    def rtype(self) -> chem.ResidueType:
        """Getter for the Residue()'s chem.ResidueType value."""
        return self.rtype_

    def get_name(self) -> str:
        """Getter for the Residue()'s name."""
        return self.name

    # === Getter-Prop (cpy/new) ===
    def num_atoms(self) -> int:
        """Number of atoms in the Residue."""
        return len(self.atoms_)

    def clone(self) -> Residue:
        """Creates a deepcopy of self."""
        return deepcopy(self)

    # === Checker === 
    def is_empty_chain(self) -> bool:
        """Whether the Residue() lacks a chain identitifer. NOT a chain object accessor."""
        return not len(self.chain_.strip())

    def is_ligand(self) -> bool:
        """Checks if Residue() is a ligand. Inherited by children."""
        return False

    def is_metal(self) -> bool:
        """Checks if Residue() is a metal. Inherited by children."""
        return False

    def is_canonical(self) -> bool:
        """Checks if Residue() is canonical. Inherited by children."""
        return self.name in chem.THREE_LETTER_AA_MAPPER #@shaoqz: @imp2 why dont we have anything related to canonical in the name

    def is_rd_solvent(self) -> bool: #@shaoqz: @imp2 move to the IO class. This is PDB related.
        """Checks if the Residue() is an rd_sovlent as defined by enzy_htp.chemical.solvent.RD_SOLVENT"""
        return self.name in chem.RD_SOLVENT_LIST

    def is_rd_non_ligand(self) -> bool:
        """Checks if the Residue() is an rd_non_ligand as defined by enzy_htp.chemical.solvent.RD_NON_LIGAND_LIST"""
        return self.name in chem.RD_NON_LIGAND_LIST

    def is_sequence_equivalent(self, other: Residue) -> bool: #@shaoqz: @imp2 why do we need to compare this? Because in the comparsion of residues there often some information exist and want to compare the rest.
        """Comparator that checks for sequence same-ness."""
        return self.residue_key == other.residue_key

    # === Editor ===
    def set_chain(self, val: str) -> None:
        """Sets chain and all child atoms to new value. Re-generates the self.residue_key attribute."""
        chain = val #@shaoqz: Refine this line
        for idx in range(len(self.atoms_)):
            self.atoms_[idx].chain_id = val
        self.chain_ = val
        self.residue_key = f"{self.chain_}.{self.name}.{self.num_}"

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
        self.atoms_ = sorted(self.atoms_, key=lambda aa: aa.atom_number) #@shaoqz: maybe use .sort to keep the reference.

        for idx, aa in enumerate(self.atoms_):
            self.atoms_[idx].atom_number = idx + start #@shaoqz: why dont use aa?
        return idx + start

    def set_rtype(self, new_rtype: chem.ResidueType) -> None:
        """Setter for the Residue()'s chem.ResidueType value."""
        self.rtype_ = new_rtype

    def sort_key(self) -> Tuple[str, int]:
        """Generates the sorting key for the Residue(). Specifically, a tuple of [chain_id, res_num]."""
        return (self.chain_, self.num_)

    # === Special ===
    def __str__(self) -> str:
        """String representation that just shows residue key in format "chain_id.residue_name.residue_num" """
        return self.residue_key

    def __repr__(self) -> str:
        """String representationt that just shows residue key in format "chain_id.residue_name.residue_num" """
        return str(self)

    #region === TODO/TOMOVE ===
    def has_alt_loc(self) -> bool:
        """Checks if any of the child Atom() objects have non-emptry alt_loc values."""
        result: bool = False
        for aa in self.atoms_:
            result |= aa.has_alt_loc() #@shaoqz: cool!
        return result

    def resolve_alt_loc(self, keep: str = "first") -> None: #@shaoqz: nice we have this now!
        """Resolves the alt_loc atoms owned by the current Residue().
        Optional argument "keep" specifies the resolution method. Default
        is "first" which forces Residue() to keep lowest location. Otherwise,
        keeps the specific alt_loc identifier, assuming it exists.
        """
        if keep == "all":
            return
        mapper: Dict[str, List[Atom]] = defaultdict(list)
        for aa in self.atoms_:
            mapper[aa.alt_loc.strip()].append(aa)

        if len(mapper.keys()) == 1:
            return

        self.atoms_ = mapper[""]
        del mapper[""]

        if keep == "first":
            keys: List[str] = sorted(list(mapper.keys()))
            self.atoms_.extend(mapper[keys[0]])
        else:
            assert keep in mapper
            self.atoms_.extend(mapper[keep])

        self.atoms_.sort(key=lambda aa: aa.atom_number)

    def remove_alt_loc(self) -> None:
        """Makes the alt_loc for each of the child Atom() objects blank."""
        for idx in range(len(self.atoms_)):
            self.atoms_[idx].alt_loc = ""

    def remove_occupancy(self) -> None: #@shaoqz: @imp2 really curious about the reason of putting all fundmental methods that edit data in residue instead of atom, which the later is more intuative to me.
        """Makes the occupancy for each of the child Atom() objects blank."""
        for idx in range(len(self.atoms_)):
            self.atoms_[idx].occupancy = ""

    def min_line(self) -> int: #@shaoqz: @imp2 move to the pdb IO class
        """The lowest one-indexed line one of the Residue()'s atoms come from."""
        return self.min_line_

    def max_line(self) -> int:
        """The highest one-indexed line one of the Residue()'s atoms come from."""
        return self.max_line_

    def line_range(self) -> Tuple[int, int]:
        """A tuple of the he lowest and highest one-indexed lines that the Residue()'s atoms come from."""
        return (self.min_line(), self.max_line())

    def neighbors(self, other: Residue) -> bool: #@shaoqz: @imp2 move to IO class / change to is_neighbors
        """Checks if two residues from the same PDB are neighbors."""
        return (abs(self.min_line() - other.max_line()) == 1) or (
            abs(self.max_line() - other.min_line()) == 1
        )

    def get_pdb_lines(self, terminal: bool = False) -> List[str]:
        """Method that gets the PDB lines for all of the Atom() objects in the Residue() object."""
        result = list()
        for aa in self.atoms_:
            result.append(aa.to_pdb_line())
        if terminal:
            num = int(result[-1][6:11])
            ter_line = result[-1]
            ter_line = ter_line[:26] + " " * (len(ter_line) - 26)
            ter_line = "TER   " + ter_line[6:]
            print(result[-1][6:12])
            num = int(result[-1][6:12])
            ter_line = ter_line[0:6] + f"{num+1: >5d}     " + ter_line[16:]
            result.append(ter_line)
        return result

    #endregion
