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

    #endregion

    # ============= TODO below =================
    #region === Getters === (Attributes - accessing Structure data -  references)
    def get_atom(self) -> Atom:
        """TODO do we really need this? Providing access of deeper layer requires a key to select
        maybe just use python objects to access is a better idea.
        And these APIs are just for developers, users will have selector in the future to do selection"""
        pass

    #endregion

    #region === Getter === (Properities - derived data; wont affect Structure data - copy)
    @property
    def residue_state(self) -> List[Tuple[str, str, int]]:  #@shaoqz: @residue_key
        """Generates a list of tuples of all residues in the Structure. Format for each tuple is (one_letter_res_name, chain_id, res_index).
        This method is designed for debuging purpose"""
        result = list()
        for cname, chain in self.chain_mapper.items():
            for residue in chain.residues:
                (chain, res_name, index) = residue.residue_key.split(".")
                if residue.is_canonical():
                    result.append((chain, convert_to_one_letter(res_name), int(index)))
                elif residue.is_metal():  #@shaoqz: @imp2 any non-canonical should be using 3-letter name
                    result.append((chain, res_name, int(index)))
        return result

    @property
    def residue_keys(self) -> List[str]:
        """Generates a list of strings containing all residue_key values for all child Residue()"s"""
        result = list()
        for chain in self.chains:
            for res in chain.residues:
                result.append(res.residue_key)
        return result

    @property  # @shaoqz: @imp2 do we really need this here? residue key is nessessary since we are operating acrossing levels
    def num_residues(self) -> int:
        """Returns the number of Residue() objects contained within the current Structure()."""
        total: int = 0
        ch: Chain
        for ch in self.chains:
            total += ch.num_residues
        return total

    @property
    def chain_names(self) -> List[str]:
        """Returns a list of all the chain names for the Structure()"""
        return list(self.chain_mapper.keys())

    #endregion

    #region === Checker ===
    def has_chain(self, chain_name: str) -> bool:
        """Checks if the Structure() has a chain with the specified chain_name."""
        return chain_name in self.chain_mapper

    #endregion


    #region === Editor ===
    def add_chain(self, new_chain: Chain, overwrite: bool = False) -> None:  #TODO add logic for overwriting
        """Method that inserts a new chain and then sorts the chains based on name.
        Will overwrite if Chain() with existing name already in object. #@shaoqz: add + sort = insert
        """
        new_chain_name: str = new_chain.name()
        if new_chain_name in self.chain_mapper:
            self.remove_chain(
                new_chain_name
            )  #@shaoqz: give a warning. @imp should not overwrite. Since want to add a chain to a structure with a same-naming chain is very common. (like I want to merge 2 single chain object to a dimer) A better default strategy is to insert after the chain with the same name and also record a index map.

        self.chains.append(new_chain)
        self.chain_mapper[new_chain.name()] = new_chain
        self.chains.sort(key=lambda c: c.name())

    def remove_chain(self, chain_name: str) -> None:
        """Given a chain name, removes the Chain() object form both self.chains_ and self.chain_mapper."""
        del self.chain_mapper[chain_name]
        to_remove = -1
        for idx, chain in enumerate(self.chains):
            if chain.name() == chain_name:
                to_remove = idx
                break

        if to_remove != -1:
            del self.chains[to_remove]

    def add_residue(self, new_res: Residue, chain_name:str=None) -> None:
        """Inserts a new Residue() object into the Structure(). If the exact Residue (chain_id, name, residue_id) already
        exists, the new Residue overwrites it. If the new Residue specifies a Chain() that does not exist, a new chain is made.
        """
        if chain_name is None:
            chain_name: str = new_res.chain()

        if not self.has_chain(chain_name):
            new_chain: Chain = Chain(chain_name, [new_res])  #@shaoqz: give a warning
            self.chains.append(new_chain)
            self.chain_mapper[chain_name] = new_chain
        else:
            self.chain_mapper[chain_name].add_residue(new_res)

        self.chains = list(self.chain_mapper.values())  #@shaoqz: why need this
        self.chains.sort(key=lambda c: c.name)  #@shaoqz: should this be in the 1st if block?

    def remove_residue(self, target_key: str) -> None:
        """Given a target_key str of the Residue() residue_key ( "chain_id.residue_name.residue_number" ) format,
        the Residue() is removed if it currently exists in one of the child Chain()"s. If the Chain() is empty after this
        removal, the chain is deleted."""
        (chain_name, _, _) = target_key.split(".")  #@shaoqz: why not use this in get lolll
        if self.has_chain(chain_name):
            self.chain_mapper[chain_name].remove_residue(target_key)
            if self.chain_mapper[chain_name].empty():
                self.remove_chain(chain_name)

    #endregion

    #region === Special ===
    def __bool__(self) -> bool:
        """Enables running assert Structure(). Checks if there is anything in the structure."""
        return bool(len(self.chains))

    def __eq__(self, other: Structure) -> bool:  # TODO
        """Comparison operator for other Structure() objects. Checks first if both have same chain names and then if each named chain is identical."""
        if set(self.chain_mapper.keys()) != set(other.chain_mapper.keys()):
            return False

        chain_name: Chain
        other_chain: Chain
        for chain_name, self_chain in self.chain_mapper.items():
            other_chain = other.chain_mapper[chain_name]
            if not self_chain.is_same_sequence(other_chain):
                return False
        return True  #@shaoqz: so this comparsion is only in sequence level. This does not really make sense. Different levels of comparsion
        #         is needed for a pair structure
        #         for example:
        #            - if the *coordinate* of every atom is the same
        #            - if the *topology* is the same. (presence of atoms and their connectivity)
        #            - if the *sequence* is the same. (this is basic works as the topology one as topology inside and between
        #              residues are relatively conserved, but there are exceptions, e.g. different ligand both name LIG)

    def __ne__(self, other: Structure) -> bool:
        """Negation operator for other Structure() objects. Inverstion of Structure.__eq__()."""
        return not (self == other)

    #endregion

    #region (TODO+OLD)
    # === TODO ===
    def get_connectivty_table(  #@shaoqz: ok seems not using
            self, ff="GAUSSIAN", metal_fix=1, ligand_fix=1, prepi_path=None):
        """
        get connectivity table with atom index based on "ff" settings:
        ff = GAUSSIAN  -- continuous atom index start from 1, do not seperate by chain
        -------------------
        TREATMENT
        Use original atom.id.
            chain: based on connectivity map of each atom in each residue
            metalatom:  fix1: treat as isolated atom
                        fix2: connect to donor atom (MCPB?)
            ligand: fix1: use antechamber generated prepin file to get connectivity.
        Use 1.0 for all connection.
            Amber force field in gaussian do not account in bond order. (Only UFF does.)
            Note that bond order less than 0.1 do not count in MM but only in opt redundant coordinate.
        """
        connectivty_table = ""
        # get connect for every atom in stru
        self.init_connect(metal_fix, ligand_fix, prepi_path)

        # write str in order
        # Note: Only write the connected atom with larger id
        a_id = 0
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    a_id += 1
                    cnt_line = " " + str(atom.id)
                    # san check
                    if atom.id != a_id:
                        raise Exception("atom id error.")
                    for cnt_atom in atom.connect:
                        if cnt_atom.id > atom.id:
                            cnt_line += " " + str(cnt_atom.id) + " " + "1.0"
                    connectivty_table += cnt_line + line_feed

        for lig in self.ligands:
            for atom in lig:
                a_id += 1
                cnt_line = " " + str(atom.id)
                # san check
                if atom.id != a_id:
                    raise Exception("atom id error.")
                for cnt_atom in atom.connect:
                    if cnt_atom.id > atom.id:
                        cnt_line += " " + str(cnt_atom.id) + " " + "1.0"
                connectivty_table += cnt_line + line_feed

        for atom in self.metalatoms:
            a_id += 1
            cnt_line = " " + str(atom.id)
            # san check
            if atom.id != a_id:
                raise Exception("atom id error.")
            for cnt_atom in atom.connect:
                if cnt_atom.id > atom.id:
                    cnt_line += " " + str(cnt_atom.id) + " " + "1.0"
            connectivty_table += cnt_line + line_feed

        for sol in self.solvents:
            for atom in sol:
                a_id += 1
                cnt_line = " " + str(atom.id)
                # san check
                if atom.id != a_id:
                    raise Exception("atom id error.")
                for cnt_atom in atom.connect:
                    if cnt_atom.id > atom.id:
                        cnt_line += " " + str(cnt_atom.id) + " " + "1.0"
                connectivty_table += cnt_line + line_feed

        return connectivty_table

    # === TO BE MOVE ===
    def build_ligands(self,
                      out_dir: str,
                      unique: bool = False) -> List[str]:  # TODO(qz): change reference of this this to get_file_str(stru.ligands[i])
        """Exports all the Ligand() objects in the Structure() to .pdb files.

        Args:
                out_dir: The base directory to save the .pdb files to.
                unique: Whether or not the saved .pdb files should be unique.

        Returns:
                A list of str() with paths to the exported ligand .pdb files.
        """
        result: List[str] = []
        existing: List[str] = []
        ligands: List[Ligand] = self.ligands

        for lidx, lig in enumerate(ligands):
            # TODO(CJ): add some kind of formatting for lidx
            lig_name: str = lig.name()
            out_pdb: str = f"{out_dir}/ligand_{lig_name}_{lidx}.pdb"

            if unique and lig_name in existing:
                continue

            lig.build(out_pdb)  #@shaoqz: use IO interface with different format instead
            result.append(out_pdb)
            existing.append(lig_name)

        return result

    def build_protein(
        self,
        dir,
        ft="PDB"
    ):  #@shaoqz: maybe unify these to a build sele option like pymol did? Support a grammer to indicate what should be contained in each file. But having these presets are also good.
        """
        build only protein and output under the dir
        -------
        dir: out put dir ($dir/protein.pdb)
        ft : file type / now support: PDB(default)
        """
        # make path
        if dir[-1] == "/":
            dir = dir[:-1]
        out_path = dir + "/protein.pdb"

        # write
        if ft == "PDB":
            with open(out_path, "w") as of:  #@shaoqz: same as build ligand use IO interface class
                a_id = 0
                r_id = 0
                for chain in self.chains:
                    # write chain
                    for resi in chain:
                        r_id = r_id + 1
                        for atom in resi:
                            a_id = a_id + 1  # current line index
                            line = atom.build(a_id=a_id, r_id=r_id)
                            of.write(line)
                    # write TER after each chain
                    of.write("TER" + line_feed)  #@shaoqz: @imp2 how do you solve these
                of.write("END" + line_feed)
        else:
            raise Exception("Support only PDB output now.")

        return out_path

    def build_metalcenters(self, dir, ft="PDB"):
        """
        build metalcenters only. Use for MCPB parameterization. Deal with donor residue with different protonation states.        ----------
        TODO
        """
        out_paths = []
        return out_paths

    # def get_atom_id(self): #@shaoqz: @nu
    #     """
    #     return a list of id of all atoms in the structure
    #     """
    #     atom_id_list = []
    #     for chain in self.chains_:
    #         for res in chain:
    #             for atom in res:
    #                 if atom.id == None:
    #                     raise Exception(
    #                         "Detected None in chain "
    #                         + str(chain.id)
    #                         + str(res.id)
    #                         + " "
    #                         + atom.name
    #                     )
    #                 atom_id_list.append(atom.id)
    #     for metal in self.metalatoms:
    #         if metal.id == None:
    #             raise Exception("Detected None in metal " + metal.name)
    #         atom_id_list.append(metal.id)
    #     for lig in self.ligands:
    #         for atom in lig:
    #             if atom.id == None:
    #                 raise Exception("Detected None in ligands", res.id, atom.name)
    #             atom_id_list.append(atom.id)
    #     for sol in self.solvents:
    #         for atom in sol:
    #             if atom.id == None:
    #                 raise Exception("Detected None in solvent", res.id, atom.name)
    #             atom_id_list.append(atom.id)
    #     return atom_id_list

    # def get_atom_type(self, prmtop_path): #@shaoqz: @nu
    #     """
    #     requires generate the stru using !SAME! PDB as one that generate the prmtop.
    #     """
    #     # get type list
    #     with open(prmtop_path) as f:
    #         type_list = []
    #         line_index = 0
    #         for line in f:
    #             line_index = line_index + 1  # current line
    #             if line.strip() == r"%FLAG POINTERS":
    #                 format_flag = line_index
    #             if line.strip() == r"%FLAG AMBER_ATOM_TYPE":
    #                 type_flag = line_index
    #             if "format_flag" in dir():
    #                 if line_index == format_flag + 2:
    #                     N_atom = int(line.split()[0])
    #                     del format_flag
    #             if "type_flag" in dir():
    #                 if (
    #                     line_index >= type_flag + 2
    #                     and line_index <= type_flag + 1 + ceil(N_atom / 5)
    #                 ):
    #                     for i in line.strip().split():
    #                         type_list.append(i)
    #     # assign type to atom
    #     for chain in self.chains_:
    #         for res in chain:
    #             for atom in res:
    #                 atom.type = type_list[atom.id - 1]
    #     for atom in self.metalatoms:
    #         if atom.id == None:
    #             raise Exception("Detected None in metal " + atom.name)
    #         atom.type = type_list[atom.id - 1]
    #     for lig in self.ligands:
    #         for atom in lig:
    #             if atom.id == None:
    #                 raise Exception("Detected None in ligands", res.id, atom.name)
    #             atom.type = type_list[atom.id - 1]
    #     for sol in self.solvents:
    #         for atom in sol:
    #             if atom.id == None:
    #                 raise Exception("Detected None in solvent", res.id, atom.name)
    #             atom.type = type_list[atom.id - 1]

    # def get_resi_dist(self, r1, r2, method="mass_center"): #@shaoqz: @nu
    #     """
    #     r1: residue 1. (residue_obj)
    #     r2: residue 2. (residue_obj)
    #     method: The method to narrow the residue down to a point.
    #             - mass_center
    #             - ...
    #     """
    #     if method == "mass_center":
    #         p1 = r1.get_mass_center()
    #         p2 = r2.get_mass_center()
    #     D = get_distance(p1, p2)
    #     return D
    #endregion

def compare_structures(left: Structure, right: Structure) -> Dict[str, List[str]]:
    """Compares two Structure() objects and returns a dict() of missing Residues with format:

    {"left": ["residue_key1","residue_key1",..],
     "right": ["residue_key1","residue_key1",..]
         }
    """
    result = {"left": [], "right": []}
    left_keys: Set[str] = set(left.residue_keys)
    right_keys: Set[str] = set(right.residue_keys)

    result["left"] = list(filter(lambda ll: ll not in right_keys, left_keys))
    result["right"] = list(filter(lambda rr: rr not in left_keys, right_keys))
    return result


def merge_right(left: Structure,
                right: Structure) -> Structure:  #@shaoqz: I believe there will be a bug due to the treatment of insert with same id
    """Merges Residue() and derived objects from left Structure() to right Structure(), making sure that ALL Residue() and
        Residue() derived objects from the left are in the right. Note that the reverse is not applied and that elements initially found only in right are NOT
        merged back over to left. Also not the resulting Structure() is a deepcopy and no changes are made to the original left or right objects.

    Example:
            #TODO(CJ): add this in
    """
    struct_cpy: Structure = deepcopy(right)
    # TODO(CJ): make this a method
    left.chains = sorted(left.chains, key=lambda c: c.name())
    struct_cpy.chains = sorted(struct_cpy.chains, key=lambda c: c.name())
    # this is the case where there is straight up a missing chain
    for cname, chain in left.chain_mapper.items():
        if not struct_cpy.has_chain(cname):
            struct_cpy.add_chain(deepcopy(chain))

    right_keys = struct_cpy.residue_keys  #@shaoqz: @imp2 why is this part needed as all res is stored in chain.
    for lkey in left.residue_keys:
        if lkey not in right_keys:
            struct_cpy.add_residue(left.get_residue(lkey))
    return struct_cpy
