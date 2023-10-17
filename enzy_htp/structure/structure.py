"""This class define the core data structure of EnzyHTP: Structure. Structure stands for the single enzyme structure.
As the core data class, Structure will be solely designed for **storing, accessing and editing** structural data.

For the data point of view, the enzyme structure is composed by two parts of data:
- topology (composition of atoms and their connectivity)
  there are 2 common ways to store connectivity:
    + the chain, residue division of atoms (inter-residue connectivity)
      & the canonical residue name and atom names in it (intra-residue connectivity)
      & the connectivity table for non-canonical parts
    + the connectivity table for the whole enzyme
- coordinate (atomic coordinate upon topology)

Base on this concept, Structure object holds a list of composing Chain objects which also contain Residue objects
and so that Atom objects. In each Atom object, atom name and coordinate is stored. Every level (chain, residue, atom/
coordinate/connectivity) of information can all be pulled out from the Structure object with getter methods and can be set
with setter methods. Also Structure() supports common editing methods such as add/remove children objects.

Application of Structure objects - Binding modules:
Generation:
    Note that Structure() objects SHOULD NOT be created by the user directly and instead created through different generation 
    methods from the binding StructureIO classes (e.g.: enzy_htp.structure_io.PDBParser().get_structure()) from different file 
    types and different external data structures.
Selection:
    Selection of Structural regions are handled by the Selection module.
Operation:
    Changes of the Structure data are handled by functions in the operation module. These commonly used operations of Structure 
    will than be used in scientific APIs: Preparation, Mutation, Geom Variation. And structure based descriptors are derived by 
    functions in the Energy Engine module.

### TO BE UPDATE ###
    Sub-module also contains two utility free functions for comparing and merging structural elements.
    compare_structures() and merge_right() are used for comparing and merging structural components for Structure()
    objects, respectively.

### below doc is waiting for update in the future if change ###
    NOTE: The below code assumes a trivial protein structure with two chains containing 3 and 1 residue, respectively.

    Loading a structure from PDB:
        >>> import enzy_htp
        >>> structure : enzy_htp.Structure = enzy_htp.structure_from_pdb("/path/to/pdb")

    Surveying basic information:
        >>> structure.num_chains()
        2
        >>> structure.residue_state
        [("A","ASP",1),("A","ASP",2),("A","ASP",3),("B","ASP",1)] #@shaoqz:shouldn"t it be 1-letter name?
        >>> structure.num_residues
        4

    Interfacing with Chain()'s:
        >>> structure.chains
        [<enzy_htp.structure.chain.Chain object at 0x7fdc680ac670>, 
            <enzy_htp.structure.chain.Chain object at 0x7fdc680855b0>]
        >>> structure.chain_names
        ["A", "B"]
        >>> chain_cpy : enzy_htp.Chain = structure.get_chain( "B" )
        >>> structure.remove_chain( "B" )
        >>> structure.chains
        [<enzy_htp.structure.chain.Chain object at 0x7fdc680ac670>]
        >>> structure.num_chains()
        1
        >>> structure.num_residues
        3
        >>> structure.add_chain( chain_cpy )
        >>> structure.chains
        [<enzy_htp.structure.chain.Chain object at 0x7fdc680ac670>, 
            <enzy_htp.structure.chain.Chain object at 0x7fdc680855b0>] # @shaoqz: why chain.Chain?
        >>> structure.chain_names
        ["A", "B"]
        >>> structure.num_chains()
        2	
        >>> structure.num_residues
        4	

    Interfacing with Residue()"s:
        >>> structure.residues()
        ["A.ASP.1","A.ASP.2","A.ASP.3","B.ASP.1"]   # @shaoqz: shouldn"t this be Residue objects?
        >>> structure.num_residues
        4
        >>> res_cpy : enzy_htp.Residue = structure.get_residue( "B.ASP.1" ) # @shaoqz: @imp we should not use residue name and id together as the identifier. Either id only to pinpoint or name only to batch select
        >>> structure.remove_residue( "B.ASP.1" )
        >>> structure.residues()
        ["A.ASP.1","A.ASP.2","A.ASP.3"]
        >>> structure.num_residues
        3
        >>> structure.add_residue( res_cpy )
        >>> structure.residues()
        ["A.ASP.1","A.ASP.2","A.ASP.3","B.ASP.1"]
        >>> structure.num_residues
        4

    Saving the structure:
        >>> structure.to_pdb( "/path/to/copy/of/pdb" )

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""
#TODO(CJ): add a method for changing/accessing a specific residue
from __future__ import annotations
import itertools
import os
from plum import dispatch
import string
from copy import deepcopy
import sys
from typing import List, Set, Dict, Tuple, Union
from collections import defaultdict

import enzy_htp.core.math_helper as mh
from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode
from enzy_htp.chemical import convert_to_one_letter, ResidueType

from .atom import Atom
from . import Chain
from .residue import Residue
from .metal_atom import MetalUnit
from .ligand import Ligand, residue_to_ligand
from .solvent import Solvent, residue_to_solvent
from .metal_atom import MetalUnit


class Structure(DoubleLinkedNode):
    """Protein structure.
    Designed for direct interfacing by users.
    Composed of child Chain() objects and their subsequent child Residue() objects and so on Atom() objects.
    Note: This class SHOULD NOT be created directly by users. It should be created with methods from the StructureIO module.
    Note: regardless of index assigned to residues and atoms, there is an intrinsic indexing system based on the order of 
    _children lists. This intrinsicc index can be compared with pymol's index (not id)

    Attributes:
        children/chains: List[Chain]

    Derived properties:
        residue_state : List[Tuple[str, str, int]]
        residue_keys() : List[str]
        residues : List
        metals : List[Residue]
        ligands : List[Residue]
        solvents : List
        num_chains
        num_residues
        chain_names
    """

    def __init__(self, chains: List[Chain]):
        """Constructor that takes just a list of Chain() objects as input."""
        self.set_children(chains)
        self.set_ghost_parent()
        if self.has_duplicate_chain_name():
            self.resolve_duplicated_chain_name()

    #region === Getters-attr ===
    @property
    def chains(self) -> List[Chain]:
        """alias for _children. prevent changing _children but _residues holds the same"""
        return self.get_children()

    @chains.setter
    def chains(self, val) -> None:
        """setter for chains"""
        self.set_children(val)

    @property
    def _chains(self) -> List[Chain]:
        """alias for _children. prevent changing _children but _chains holds the same"""
        return self._children

    @_chains.setter
    def _chains(self, val) -> None:
        """setter for _chains"""
        self.set_children(val)

    @property
    def chain_mapper(self) -> Dict[str, Chain]:
        mapper = {}
        if self.has_duplicate_chain_name():
            self.resolve_duplicated_chain_name()
        ch: Chain
        for ch in self._chains:
            mapper[ch.name] = ch
        return mapper

    def get_chain(self, chain_name: str) -> Union[Chain, None]:
        """Gets a chain of the given name. Returns None if the Chain() is not present."""
        return self.chain_mapper.get(chain_name, None)

    #endregion

    #region === Getter-Prop ===
    @property
    def num_chains(self) -> int:
        """Returns the number of Chain() objects in the current Structure()."""
        return len(self._chains)

    @property
    def chain_names(self) -> List[str]:
        """Returns a list of chain names"""
        return list(map(lambda x: x.name, self._chains))

    def _legal_new_chain_names(self) -> List[str]:
        """
        Small helper method that determines the legal chain names that are not the same as existing ones.
        Uses all of the available 26 capitalized letters and returns in reverse order. If all characters are used
        use numbers up to 500 as place holders. Returns an empty list when all are occupied. 
        Return:
            the legal chain name list
        """
        result = list(string.ascii_uppercase) + list(map(str, range(500)))
        result = list(filter(lambda s: s not in self.chain_names, result))
        return result

    @property
    def num_residues(self) -> int:
        """Returns the number of Residue() objects in the current Structure()."""
        return len(self.residues)

    @property
    def residue_indexes(self) -> int:
        """Returns the number of Residue() objects in the current Structure()."""
        return list(map(lambda x: x.idx, self.residues))

    @property
    def residues(self) -> List[Residue]:
        """Return a list of the residues in the Structure() object sorted by (chain_id, residue_id)"""
        result = list(itertools.chain.from_iterable(self._chains))
        result.sort(key=lambda r: r.key())
        return result

    @property
    def residue_mapper(self) -> Dict[Tuple[str, int], Residue]:
        """return a mapper of {(chain_id, residue_idx): Residue (reference)}"""
        result = {}
        for residue in self.residues:
            result[residue.key()] = residue
        return result

    def find_residue_name(self, name) -> List[Residue]:
        """find residues base on its name. Return a list of matching residues"""
        result = list(filter(lambda r: r.name == name, self.residues))
        return result

    def find_residue_with_key(self, key: Tuple[str, int]) -> Union[Residue, None]:
        """find residues base on its (chain_id, idx). Return the matching residues"""
        result = list(filter(lambda r: r.key() == key, self.residues))
        if len(result) == 0:
            _LOGGER.info(f"Didn't find any residue with key: {key}")
            return None
        if len(result) > 1:
            _LOGGER.warning(f"More than 1 residue with key: {key}. Only the first one is used. Check those residues:")
            for i in result:
                _LOGGER.warning(f"    {i}")
        return result[0]

    @property
    def atoms(self) -> List[Atom]:
        """Accessor to get the atoms in the enzyme as a list of Atom objects."""
        result = []
        for chain in self.chains:
            for residue in chain:
                result.extend(residue.atoms)
        return result

    def find_atoms_in_range(self, center: Union[Atom, Tuple[float, float, float]], range_distance: float) -> List[Atom]:
        """find atoms in {range} of {center}. return a list of atoms found"""
        result = []
        for atom in self.atoms:
            if atom.distance_to(center) <= range_distance:
                result.append(atom)
        return result

    def find_idx_atom(self, atom_idx: int) -> Atom:
        """find atom base on its idx. return a reference of the atom."""
        result = list(filter(lambda a: a.idx == atom_idx, self.atoms))
        if not result:
            _LOGGER.info(f"found 0 atom with index: {atom_idx}")
        if len(result) > 1:
            _LOGGER.warning(f"found {len(result)} atoms with index: {atom_idx}! only the 1st one is used. consider sort_everything()")
        return result[0]

    def find_idxes_atom_list(self, atom_idx_list: int) -> List[Atom]:
        """find atom base on its idx. return a list reference of the atoms."""
        result = []
        for idx in atom_idx_list:
            result.append(self.find_idx_atom(idx))
        return result

    @property
    def ligands(self) -> List[Ligand]:
        """Filters out the ligand Residue()"s from the chains in the Structure()."""
        result: List[Ligand] = list()
        for chain in self.chains:
            result.extend(list(filter(lambda r: r.is_ligand(), chain)))
        return result

    @property
    def solvents(self) -> List[Solvent]:
        """return all solvents hold by current Structure()"""
        result: List[Solvent] = []
        for chain in self.chains:
            result.extend(list(filter(lambda r: r.is_solvent(), chain)))
        return result

    @property
    def metals(self) -> List[MetalUnit]:
        """Filters out the metal Residue()"s from the chains in the Structure()."""
        result: List[Residue] = []
        for chain in self.chains:
            result.extend(list(filter(lambda r: r.is_metal(), chain.residues)))
        return result

    @property
    def metalcenters(self) -> List[MetalUnit]:
        """Filters out the metal coordination center Residue()"s from the chains 
        in the Structure()."""
        result: List[Residue] = []
        for chain in self.chains:
            result.extend(list(filter(lambda r: r.is_metal_center(), chain.residues)))
        return result

    @property
    def polypeptides(self) -> List[Chain]:
        """return the peptide part of current Structure() as a list of chains"""
        result: List[Chain] = list(filter(lambda c: c.is_polypeptide(), self._chains))
        return result

    @property
    def sequence(self) -> Dict[str, str]:
        """return a dictionary of {'chain_id':'chain_sequence'}"""
        result = {}
        self.sort_chains()
        ch: Chain
        for ch in self._chains:
            result[ch.name] = ch.sequence
        return result

    def init_connect(self,
                     ligand_fix: int = "antechamber",
                     metal_fix: int = "isolate",
                     ncaa_fix: int = "antechamber",
                     solvent_fix: int = "caa") -> None:
        """
        Initiate connectivity for the Structure.
        Save the connectivity to self._connect of each Atom().
        Args:
            ligand_fix:
                the method that determines connectivity for ligand. (see details below)
            metal_fix:
                the method that determines connectivity for metal. (see details below)
            ncaa_fix:
                the method that determines connectivity for ncaa. (see details below)
            solvent_fix:
                the method that determines connectivity for solvent. (see details below)

        Details:
            polypeptide:
                using documented connectivity for each canonical amino acid from Amber
                library.
            ligand:
                fix = "antechamber":
                    use antechamber to generated connectivity and read from prepin file.
                    (according to https://ambermd.org/doc/prep.html the coordniate line will
                    always start at the 11th line after 3 DUMM.)
            metalatom:
                fix = "isolate": treat as isolated atom
                TODO fix = "mcpb": connect to donor atom (MCPB?)
            non-canonical residues:
                fix = "antechamber": same as above in ligand.
            solvent:
                fix = "caa": same as polypeptide part
        """
        # TODO(qz)(high_prior) finish all the called functions
        self.init_connect_for_polypeptides(ncaa_fix=ncaa_fix)
        self.init_connect_for_ligands(method=ligand_fix)
        self.init_connect_for_metals(method=metal_fix)
        self.init_connect_for_solvents(method=solvent_fix)
        # check if all atoms are connected
        for atm in self.atoms:
            if not atm.is_connected():
                _LOGGER.error(f"Atom {atm} doesn't have connect record after initiation.")
                sys.exit(1)

    def init_connect_for_polypeptides(self, ncaa_fix: str):
        """initiate connectivity for polypeptides."""
        for chain in self.polypeptides:
            for res in chain.residues:
                if res.is_noncanonical():
                    res.init_connect_ncaa(method=ncaa_fix)
                else:
                    for atom in res.atoms:
                        atom.init_connect_in_caa()

    def init_connect_for_ligands(self, method: str):
        """"""
        # TODO(qz)(high_prior)
        # TODO generate prepi by itself and store it to a global database path so that other
        # process in the same workflow can share the generated file.
        support_method_list = ["caa"]
        # ligand
        # init
        for lig in self.ligands:
            for atom in lig:
                atom.connect = []
        # fix 1
        if ligand_fix == 1:
            for lig in self.ligands:
                # read prepin for each ligand
                with open(prepi_path[lig.name]) as f:
                    line_id = 0
                    if_loop = 0
                    for line in f:
                        line_id += 1
                        if line.strip() == "":
                            if if_loop == 1:
                                # switch off loop and break if first blank after LOOP encountered
                                if_loop = 0
                                break
                            continue
                        if if_loop:
                            lp = line.strip().split()
                            lig._find_atom_name(lp[0]).connect.append(lig._find_atom_name(lp[1]))
                            continue
                        # loop connect starts at LOOP
                        if line.strip() == "LOOP":
                            if_loop = 1
                            continue
                        # coord starts at 11th
                        if line_id >= 11:
                            lp = line.strip().split()
                            atom_id = int(lp[0]) - 3
                            atom_cnt = int(lp[4]) - 3
                            if atom_cnt != 0:
                                lig[atom_id - 1].connect.append(lig[atom_cnt - 1])
                                lig[atom_cnt - 1].connect.append(lig[atom_id - 1])
        if method not in support_method_list:
            _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

    def init_connect_for_metals(self, method: str):
        """initiate connectivity for metals in the structure"""
        for metal in self.metals:
            metal.init_connect(method)

    def init_connect_for_solvents(self, method: str):
        """initiate connectivity for solvents in the structure"""
        for sol in self.solvents:
            sol.init_connect(method)

    #endregion

    #region === Checker ===
    def has_charges(self) -> bool:
        """Checks if the current Structure has charges for all atoms.

        Returns:
            Whether the Structure object has non-None charges for all atoms.
        """
        for aa in self.atoms:
            if aa.charge is None:
                return False
        return True

    def has_duplicate_chain_name(self) -> bool:
        """check if self._chain have duplicated chain name
        give warning if do."""
        existing_c_id = []
        ch: Chain
        for ch in self._chains:
            if ch.name in existing_c_id:
                _LOGGER.warning(f"Duplicate chain names detected in Structure obj during {sys._getframe().f_back.f_code.co_name}()! ")
                return True
            existing_c_id.append(ch.name)
        return False

    def is_idx_subset(self, target_stru: Structure) -> bool:
        """
        check if residue indexes of the target_stru is a subset of self
        by subset it means:
        1. name of chains in target_stru is contain in self
        2. for each chain, residue indexes in target chain is contained in correponding
           chain from self.
        Args:
            target_stru: the target structure to compare with self
        Returns:
            Boolean reflect if the target structure is a index subset of self
        """
        trgt_ch: Chain
        for trgt_ch in target_stru:
            if trgt_ch.name not in self.chain_mapper:
                _LOGGER.info(f"current stru {list(self.chain_mapper.keys())} doesnt contain chain: {trgt_ch} from the target stru")
                return False

            self_ch = self.chain_mapper[trgt_ch.name]
            self_ch_resi_idxes = self_ch.residue_idxs
            res: Residue
            for res in trgt_ch:
                if res.idx not in self_ch_resi_idxes:
                    _LOGGER.info(f"current stru chain {self_ch} doesnt contain residue: {res} of the target stru")
                    return False
        return True

    def has_chain(self, chain_name: str) -> bool:
        """Checks if the Structure() has a chain with the specified chain_name."""
        return chain_name in self.chain_mapper
    #endregion

    #region === Editor ===
    def sort_chains(self) -> None:
        """
        sort children chains with their chain name
        sorted is always better than not but Structure() is being lazy here
        """
        self._chains.sort(key=lambda x: x.name)

    def sort_everything(self) -> None:
        """sort all object in structure"""
        self.sort_chains()
        for ch in self._chains:
            ch.sort_residues()
            for res in ch:
                res: Residue
                res.sort_atoms()

    def renumber_atoms(self, sort_first: bool = True) -> None:
        """give all atoms in current structure a new index start from 1.
        Will not give a idx_change_map since idx only matter in Residue level"""
        if sort_first:
            self.sort_everything()
        _LOGGER.info("renumbering atoms")
        a_id = 1
        for atom in self.atoms:
            if atom.idx != a_id:
                _LOGGER.debug(f"changing atom {atom.idx} -> {a_id}")
            atom.idx = a_id
            a_id += 1

    def resolve_duplicated_chain_name(self) -> None:
        """resolve for duplicated chain name in self.chains_
        A. insert the chain to be one character after the same one if they are different
           and move the rest chain accordingly
        B. give error if they are the same (in coordinate)"""
        mapper = {}
        if_rename = 0
        for ch in self.chains:
            if ch.name in mapper:
                if ch.is_same_coord(mapper[ch.name]):
                    _LOGGER.error("Duplicate chain (same coordinate) detected in Structure obj! Exiting... ")
                    sys.exit(1)
                new_name = chr(ord(ch.name) + 1)  # TODO find a way
                ch.name = new_name
                if_rename = 1
            mapper[ch.name] = ch
        if if_rename:
            _LOGGER.warning("Resolved duplicated chain (different ones) name by renaming.")

    @dispatch
    def add(self, target: Chain,
            overwrite: bool = False, sort: bool = False) -> None:
        """Method that inserts a new chain.
        Args:
            overwrite: if overwrite when chain has same name. if not overwrite the new
                chain will be assigned with a new name.
            sort: if sort after adding"""
        new_chain_name: str = target.name
        if new_chain_name in self.chain_mapper:
            if overwrite:
                self.chain_mapper[new_chain_name].delete_from_parent()
            else:
                target.name = self._legal_new_chain_names()[0]
        target.parent = self
        self.chains.append(target)
        if sort:
            self.sort_chains()

    @dispatch
    def add(self, target: List[Chain], # pylint: disable=function-redefined
            overwrite: bool = False, sort: bool = False) -> None: 
        """add a list of chains into the structure."""
        for ch in target:
            self.add(ch, overwrite, sort=False)
        if sort:
            self.sort_chains()

    @dispatch
    def add(self, target: Residue, # pylint: disable=function-redefined
            sort: bool = False, chain_name:str=None) -> None: 
        """add a residue into the structure."""
        res_type = target.rtype
        if res_type in [ResidueType.CANONICAL,
                        ResidueType.NONCANONICAL]:
            _LOGGER.error("adding amino acid into the structure needs to be chain-specific. Please use Chain().add()")
            raise NameError
        elif res_type in [ResidueType.LIGAND,
                          ResidueType.METAL,
                          ResidueType.SOLVENT,
                          ResidueType.UNKNOWN]:
            # always make a new chain as they are not covalently connected
            if not chain_name:
                chain_name = self._legal_new_chain_names()[0]
            
            new_chain = Chain(
                name=chain_name,
                residues=[target,],
                parent=self,)
            self.chains.append(new_chain)
        else:
            _LOGGER.error(f"residue type {res_type} not supported")
            raise TypeError
        if sort:
            for ch in self:
                ch.sort_residues()

    @dispatch
    def add(self, target: List[Residue], # pylint: disable=function-redefined
            sort: bool = False) -> None: 
        """add a list of residues into the structure."""
        for res in target:
            self.add(res, False)
        if sort:
            for ch in self:
                ch.sort_residues()

    #endregion

    #region === Special ===
    def __str__(self):
        """
        a string representation of Structure()
        Goal:
            show it is a Structure obj
            show the data in current instance
        """
        out_line = [
            f"<Structure object at {hex(id(self))}>",
        ]
        out_line.append("Structure(")
        out_line.append(f"chains: (sorted, original {list(self.chain_mapper.keys())})")
        for ch in sorted(self._chains, key=lambda x: x.name):
            ch: Chain
            out_line.append(f"    {ch.name}({ch.chain_type}): residue: {ch.residue_idx_interval()} atom_count: {ch.num_atoms}")
        out_line.append(")")
        return os.linesep.join(out_line)

    def __getitem__(self, key: Union[int, str]):
        """support dictionary/list-like access"""
        if isinstance(key, int):
            return super().__getitem__(key)
        if isinstance(key, str):
            return self.chain_mapper[key]
        raise KeyError("Structure() getitem only take int or str as key")

    def __delitem__(self, key: Union[int, str]):
        """support dictionary like delete"""
        if isinstance(key, int):
            super().__delitem__(key)
        if isinstance(key, str):
            self.chain_mapper[key].delete_from_parent()
        raise KeyError("Structure() delitem only take int or str as key")

    def __bool__(self) -> bool:
        """Enables running assert Structure(). Checks if there is anything in the structure."""
        return bool(len(self.chains) != 0)

    def __eq__(self, other: Structure) -> bool:
        """Structure comparsion is a multi-demension task. Vaguely asking for comparing just structures is not allowed."""
        _LOGGER.error("Vaguely asking for comparing just structures is not allowed. Please use Structure().same_xxx. (xxx stands for a specific demension)")
        raise NameError

    def __ne__(self, other: Structure) -> bool:
        """Structure comparsion is a multi-demension task. Vaguely asking for comparing just structures is not allowed."""
        _LOGGER.error("Vaguely asking for comparing just structures is not allowed. Please use Structure().same_xxx. (xxx stands for a specific demension)")
        raise NameError
    #endregion

    @dispatch
    def _(self):
        """
        dummy method for dispatch
        """
        pass
