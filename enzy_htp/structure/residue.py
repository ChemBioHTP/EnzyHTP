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

from enzy_htp.core.doubly_linked_tree import DoubleLinkNode

from .atom import Atom
from enzy_htp.core import _LOGGER
import enzy_htp.chemical as chem


class Residue(DoubleLinkNode):
    """Most common functional unit in enzy_htp. Made up of Atom() objects and can be either
        canonical, non-canonical, solvent or a metal center.

    Attributes:
        (nessessary)
        _children/_atoms : A list of Atom() objects that make up the Residue().
        _name : Residue name.
        _idx : The index of the Residue within the chain.
        _parent/chain : Parent chain name.
        _rtype : The ResidueType of the Residue().
    """

    def __init__(self, residue_idx: int, residue_name: str, atoms: List[Atom], parent=None):
        """Constructor for the Residue() object"""
        self._name = residue_name
        self._idx = residue_idx
        self._rtype = chem.ResidueType.UNKNOWN
        self.set_parent(parent)
        self.set_children(atoms)

        self._atoms = self._children # add an alias
        
    #region === Getter-Attr (ref) ===
    @property
    def atoms(self) -> List[Atom]:
        """Returns a list of all Atom() objects that the Residue() "owns" """
        return self.get_children()
    @atoms.setter
    def atoms(self, val):
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
        """Getter for the Residue()'s name."""
        return self._name
    @name.setter
    def name(self, val):
        self._name = val
    #endregion
    
    #region === Getter-Prop (cpy/new) ===
    @property
    def key(self, if_name: bool = False):
        '''
        an unique tuple to indicate a specific residue in the enzyme
        '''
        if if_name:
            return (self.chain.id, self._idx, self._name)
        return (self.chain.id, self._idx)

    def num_atoms(self) -> int:
        """Number of atoms in the Residue."""
        return len(self._atoms)

    def clone(self) -> Residue:
        """Creates a deepcopy of self."""
        return deepcopy(self)
    #endregion

    #region === Checker === 
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
        return self.name in chem.THREE_LETTER_AA_MAPPER

    def is_solvent(self) -> bool: #@shaoqz: @imp2 move to the IO class. This is PDB related.
        """Checks if the Residue() is an rd_sovlent as defined by enzy_htp.chemical.solvent.RD_SOLVENT"""
        return self.name in chem.RD_SOLVENT_LIST

    def is_rd_non_ligand(self) -> bool:
        """Checks if the Residue() is an rd_non_ligand as defined by enzy_htp.chemical.solvent.RD_NON_LIGAND_LIST"""
        return self.name in chem.RD_NON_LIGAND_LIST

    def is_sequence_eq(self, other: Residue) -> bool: #@shaoqz: @imp2 why do we need to compare this? Because in the comparsion of residues there often some information exist and want to compare the rest.
        """Comparator that checks for sequence same-ness."""
        return self.key == other.key
    #endregion

    # === Editor ===
    def set_chain(self, val: str) -> None:
        """Sets chain and all child atoms to new value. Re-generates the self.residue_key attribute."""
        chain = val #@shaoqz: Refine this line
        for idx in range(len(self._atoms)):
            self._atoms[idx].chain_id = val
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
        self._atoms = sorted(self._atoms, key=lambda aa: aa.atom_number) #@shaoqz: maybe use .sort to keep the reference.

        for idx, aa in enumerate(self._atoms):
            self._atoms[idx].atom_number = idx + start #@shaoqz: why dont use aa?
        return idx + start

    def set_rtype(self, new_rtype: chem.ResidueType) -> None:
        """Setter for the Residue()'s chem.ResidueType value."""
        self._rtype = new_rtype

    def sort_key(self) -> Tuple[str, int]:
        """Generates the sorting key for the Residue(). Specifically, a tuple of [chain_id, res_num]."""
        return (self.chain_, self.num_)

    # === Special ===
    def __str__(self) -> str:
        """String representation that just shows residue key in format "chain_id.residue_name.residue_num" """
        return f'Residue({self.idx}, {self.name}, Atoms:{len(self.atoms)}, {self._parent})'

    def __repr__(self) -> str:
        """String representationt that just shows residue key in format "chain_id.residue_name.residue_num" """
        return str(self)

    #region === TODO/TOMOVE ===
    def get_pdb_lines(self, terminal: bool = False) -> List[str]:
        """Method that gets the PDB lines for all of the Atom() objects in the Residue() object."""
        result = list()
        for aa in self._atoms:
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
