"""The module define the class of selection of the part of a Structure
The instance normally should be made by enzy_htp.structure.structure_selection.general.py::select_stru

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-15
"""

import os
from typing import List, Union

from enzy_htp.core import _LOGGER
from enzy_htp.core.general import get_interval_str_from_list
from enzy_htp.core.exception import WrongTopology
from enzy_htp.structure import Structure, Atom, Residue, Chain


class StruSelection:
    """Structure selection result.
    Designed for access of the selection result because the result normally wont be a 
    simple reference of e.g. a Chain of the Structure but a composition of references.

    Attributes:
        atoms: a list of Atom() which are references from the Structure of selection

    Derived properties:
        involved_residues
        involved_chains
        atom_idx
        involved_residue_idx
        involved_chain_name"""

    def __init__(self, atoms: List[Atom]) -> None:
        """carry a list of Atom references of the Structure of selection"""
        self._atoms = atoms

    # region == getter attr ==
    @property
    def atoms(self) -> List[Atom]:
        """getter for _atoms."""
        return self._atoms

    @atoms.setter
    def atoms(self, val):
        """Normally this should not be set"""
        _LOGGER.error("StruSelection atoms should not change once constructed")
        raise Exception("StruSelection atoms should not change once constructed")
    # endregion

    # region == getter prop ==
    @property
    def involved_residues(self) -> List[Residue]:
        """return a list of residues involved in holding atoms"""
        return list(set(map(lambda x: x.residue, self.atoms)))

    @property
    def involved_chains(self) -> List[Chain]:
        """return a list of chain involved in holding atoms"""
        return list(set(map(lambda x: x.chain, self.involved_residues)))

    def atom_idx(self, if_interval_str: bool = False) -> Union[List[int], str]:
        """return a list of atom idx of holding"""
        result = list(map(lambda x: x.idx, self.atoms))
        if if_interval_str:
            result = get_interval_str_from_list(result)
        return result

    def involved_residue_idx(self, if_interval_str: bool = False) -> Union[List[int], str]:
        """return a list of residue idx involved in holding atoms"""
        result = list(map(lambda x: x.idx, self.involved_residues))
        if if_interval_str:
            result = get_interval_str_from_list(result)
        return result

    def involved_chain_name(self) -> List[str]:
        """return a list of chain name involved in holding atoms"""
        result = list(map(lambda x: x.name, self.involved_chains))
        return result
    # endregion

    # region == checker ==
    def check_consistent_topology(self) -> None:
        """check whether Atom()s are from the same Structure()
        raise an error if not.
        TODO this might suggestion StruSelection() and the attribute of
        StructureConstraint should have a common base class?"""
        top = self.atoms[0].root()
        for atom in self.atoms:
            current_top = atom.root()
            if not isinstance(current_top, Structure):
                _LOGGER.error(
                    f"Topology should be Structure(). ({atom} has {current_top})")
                raise WrongTopology

            if current_top is not top:
                _LOGGER.error(
                    "Atom()s in StructureConstraint() "
                    "have inconsistent topology! "
                    f"({self.atoms[0]} has {top} --vs-- {atom} has {current_top})")
                raise WrongTopology
    # endregion

    # region == special ==
    def __str__(self) -> str:
        """the string represetation of the selection"""
        result_line = []
        result_line.append(f"In atoms: {self.atom_idx(if_interval_str=True)}")
        result_line.append(f"Involve residues: {self.involved_residue_idx(if_interval_str=True)}")
        result_line.append(f"Involve chains: {','.join(self.involved_chain_name())}")

        return f"{os.linesep}".join(result_line)
    # endregion

    # == constructor ==
    @classmethod
    def from_atom_idx_list(cls, stru: Structure, atom_idx_list: List[int]):
        """create a StruSelection() from a list of atom indexes"""
        atom_obj_refs = stru.find_idxes_atom_list(atom_idx_list)
        return cls(atom_obj_refs)
