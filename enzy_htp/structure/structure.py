"""Definition for the Structure class. Structure objects represent a single protein/enzyme system. They
are composed of Chain objects and their respecitve residues. The majority of EnzyHTP's structural 
manipulations are carried out through this class. Note that Structure() objects SHOULD NOT be created
by the user directly and instead creted through enzy_htp.structure_from_pdb(). Common Structure() operations 
include the insertion/deletion of child Residue() and Chain() objects. 

Sub-module also contains two utility free functions for comparing and merging structural elements.
compare_structures() and merge_right() are used for comparing and merging structural components for Structure()
objects, respectively.

NOTE: The below code assumes a trivial protein structure with two chains containing 3 and 1 residue, respectively.

Loading a structure from PDB:
	>>> import enzy_htp
	>>> structure : enzy_htp.Structure = enzy_htp.structure_from_pdb("/path/to/pdb")

Surveying basic information:
    >>> structure.num_chains()
	2
	>>> structure.residue_state()
	[('A','ASP',1),('A','ASP',2),('A','ASP',3),('B','ASP',1)] #@shaoqz:shouldn't it be 1-letter name?
	>>> structure.num_residues()
	4

Interfacing with Chain()'s:
    >>> structure.chains()
    [<enzy_htp.structure.chain.Chain object at 0x7fdc680ac670>, 
	    <enzy_htp.structure.chain.Chain object at 0x7fdc680855b0>]
    >>> structure.chain_names()
    ['A', 'B']
    >>> chain_cpy : enzy_htp.Chain = structure.get_chain( 'B' )
    >>> structure.remove_chain( 'B' )
    >>> structure.chains()
    [<enzy_htp.structure.chain.Chain object at 0x7fdc680ac670>]
	>>> structure.num_chains()
	1
	>>> structure.num_residues()
	3
	>>> structure.insert_chain( chain_cpy )
    >>> structure.chains()
    [<enzy_htp.structure.chain.Chain object at 0x7fdc680ac670>, 
	    <enzy_htp.structure.chain.Chain object at 0x7fdc680855b0>] # @shaoqz: why chain.Chain?
    >>> structure.chain_names()
    ['A', 'B']
	>>> structure.num_chains()
    2	
	>>> structure.num_residues()
    4	

Interfacing with Residue()'s:
    >>> structure.residues()
    ['A.ASP.1','A.ASP.2','A.ASP.3','B.ASP.1']   # @shaoqz: shouldn't this be Residue objects?
	>>> structure.num_residues()
	4
	>>> res_cpy : enzy_htp.Residue = structure.get_residue( 'B.ASP.1' ) # @shaoqz: @imp we should not use residue name and id together as the identifier. Either id only to pinpoint or name only to batch select
    >>> structure.remove_residue( 'B.ASP.1' )
    >>> structure.residues()
    ['A.ASP.1','A.ASP.2','A.ASP.3']
	>>> structure.num_residues()
	3
	>>> structure.insert_residue( res_cpy )
    >>> structure.residues()
    ['A.ASP.1','A.ASP.2','A.ASP.3','B.ASP.1']
	>>> structure.num_residues()
	4

Saving the structure:
    >>> structure.to_pdb( "/path/to/copy/of/pdb" )


Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""
#TODO(CJ): add a method for changing/accessing a specific residue
from __future__ import annotations

import string
import pandas as pd
from copy import deepcopy
from typing import List, Set, Dict, Tuple
from collections import defaultdict
from biopandas.pdb import PandasPdb

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.chemical import one_letters_except, convert_to_one_letter

from .atom import Atom
from . import Chain
from .residue import Residue
from .metal_atom import MetalAtom
from .ligand import Ligand, residue_to_ligand
from .solvent import Solvent, residue_to_solvent
from .metal_atom import MetalAtom


class Structure:
    """High level representation of protein/enzyme structure designed for direct interfacing by users.
    Composed of child Chain() objects and their subsequent child Residue() objects. High level wrappers
    in the Structure() class enable direct manipulation of both Chain() and Residue() state while abstracting
    these details away from the user.
    Note: This class SHOULD NOT be created directly by users. It should be created with the enzy_htp.structure.structure_from_pdb() method.

    Attributes:
            chains_ : A list of Chain objects contained within the structure.
            chain_mapper : A dict() which maps chain id to the chain object.
    """

    def __init__(self, chains: List[Chain]):
        """Constructor that takes just a list of Chain() objects as input."""
        self.chains_ = chains #@shaoqz: why tailing _ #@shaoqz: @imp2 why do we need this instead of just a mapper and the getter of chains
        self.chain_mapper = dict()
        ch: Chain
        for ch in chains:
            if ch.name() in self.chain_mapper:
                _LOGGER.error(
                    "Duplicate chain names given to Structure() object. Exiting..."
                )
                exit(1)
            self.chain_mapper[ch.name()] = ch
    #@shaoqz: @prefer use property deco to distinguish more with actual non-getter methods
    def residue_state(self) -> List[Tuple[str, str, int]]: #@shaoqz: @imp2 why this name. Maybe more doc or a better name to explain the concept (where to use it) of this function. 
        """Generates a list of tuples of all residues in the Structure. Format for each tuple is (chain_id, one_letter_res_name, res_index)."""
        result = list()
        for cname, chain in self.chain_mapper.items():
            for residue in chain.residues(): #@shaoqz: make a getitem
                (chain, res_name, index) = residue.residue_key.split(".")
                if residue.is_canonical():
                    result.append((chain, convert_to_one_letter(res_name), int(index)))
                elif residue.is_metal(): #@shaoqz: any non-canonical should be using 3-letter name
                    result.append((chain, res_name, int(index)))
        return result

    def residue_keys(self) -> List[str]:
        """Generates a list of strings containing all residue_key values for all child Residue()'s"""
        result = [] #@shaoqz: unify the way to init list
        for chain in self.chains_:
            for res in chain.residues():
                result.append(res.residue_key)
        return result

    def get_metals(self) -> List[Residue]: #@shaoqz: should be MetalAtom
        """Filters out the metal Residue()'s from the chains in the Structure()."""
        result: List[Residue] = list()
        for chain in self.chains_:
            result.extend(list(filter(lambda r: r.is_metal(), chain.residues())))
        return result

    def get_ligands(self) -> List[Residue]:
        """Filters out the ligand Residue()'s from the chains in the Structure()."""
        result: List[Residue] = list()
        for chain in self.chains_:
            result.extend(list(filter(lambda r: r.is_ligand(), chain.residues())))
        return result

    def remove_chain(self, chain_name: str) -> None:
        """Given a chain name, removes the Chain() object form both self.chains_ and self.chain_mapper."""
        del self.chain_mapper[chain_name]
        to_remove = -1
        for idx, chain in enumerate(self.chains_):
            if chain.name() == chain_name:
                to_remove = idx
                break

        if to_remove != -1:
            del self.chains_[to_remove]

    def insert_chain(self, new_chain: Chain) -> None:
        """Method that inserts a new chain and then sorts the chains based on name.
        Will overwrite if Chain() with existing name already in object. #@shaoqz: add + sort = insert
        """
        new_chain_name: str = new_chain.name()
        if new_chain_name in self.chain_mapper:
            self.remove_chain(new_chain_name) #@shaoqz: give a warning. @imp should not overwrite. Since want to add a chain to a structure with a same-naming chain is very common. (like I want to merge 2 single chain object to a dimer) A better default strategy is to insert after the chain with the same name and also record a index map.

        self.chains_.append(new_chain)
        self.chain_mapper[new_chain.name()] = new_chain
        self.chains_.sort(key=lambda c: c.name())

    def num_chains(self) -> int:
        """Returns the number of Chain() objects in the current Structure()."""
        return len(self.chains_)

    def num_residues(self) -> int:
        """Returns the number of Residue() objects contained within the current Structure()."""
        total: int = 0
        ch: Chain
        for ch in self.chains_:
            total += ch.num_residues()
        return total

    def get_chain(self, chain_name: str) -> Union[Chain, None]:
        """Gets a chain of the given name. Returns None if the Chain() is not present."""
        return self.chain_mapper.get(chain_name, None)

    def has_chain(self, chain_name: str) -> bool:
        """Checks if the Structure() has a chain with the specified chain_name."""
        return chain_name in self.chain_mapper

    def residues(self) -> List[Residue]:
        """Creates a deep copied list of the residues in the Structure() object.""" #@shaoqz: is it really a deep copy? it is not supposed to be anyway. One of the main task of Structure is to edit (record the change of ) the structure as needed 
        result = list()
        for ch in self.chains_:
            result.extend(ch.residues())
        result.sort(key=lambda r: r.sort_key())
        return result

    def to_pdb(self, out_path: str) -> None: #@shaoqz: @imp2 move to the PDB interface in the future
        """Saves the structure to the specified file in the PDB file format."""
        lines = list()
        a_idx = 1 #@shaoqz: not nessessary need to number from 1 everytime. It's better to leave a few options to meet different need.
        for cname, chain in self.chain_mapper.items():
            a_idx = chain.renumber_atoms(a_idx)
            lines.extend(chain.get_pdb_lines())
            a_idx += 1
        lines.append("END")
        fs.write_lines(out_path, lines) #@shaoqz: @imp2 need to return a mapping of new indexes and old indexes

    def get_solvents(self):
        # TODO(CJ)
        return []

    def build_ligands(self, out_dir: str, unique: bool = False) -> List[str]: #@shaoqz: @imp2 add more file format option
        """Exports all the Ligand() objects in the Structure() to .pdb files.

        Args:
                out_dir: The base directory to save the .pdb files to.
                unique: Whether or not the saved .pdb files should be unique.

        Returns:
                A list of str() with paths to the exported ligand .pdb files.
        """
        result: List[str] = []
        existing: List[str] = []
        ligands: List[Ligand] = self.get_ligands()

        for lidx, lig in enumerate(ligands):
            # TODO(CJ): add some kind of formatting for lidx
            lig_name: str = lig.get_name()
            out_pdb: str = f"{out_dir}/ligand_{lig_name}_{lidx}.pdb"

            if unique and lig_name in existing:
                continue

            lig.build(out_pdb) #@shaoqz: use IO interface with different format instead
            result.append(out_pdb)
            existing.append(lig_name)

        return result

    def build_protein(self, dir, ft="PDB"): #@shaoqz: maybe unify these to a build sele option like pymol did? Support a grammer to indicate what should be contained in each file. But having these presets are also good.
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
            with open(out_path, "w") as of: #@shaoqz: same as build ligand use IO interface class
                a_id = 0
                r_id = 0
                for chain in self.chains_:
                    # write chain
                    for resi in chain:
                        r_id = r_id + 1
                        for atom in resi:
                            a_id = a_id + 1  # current line index
                            line = atom.build(a_id=a_id, r_id=r_id)
                            of.write(line)
                    # write TER after each chain
                    of.write("TER" + line_feed) #@shaoqz: @imp2 how do you solve these
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

    def get_connect(self, metal_fix=1, ligand_fix=1, prepi_path=None): #@shaoqz: @nu ### TODO rewrite doc # also change a name this is not a getter but to generate sth.
        """
        get connectivity 
        -----------------
        TREATMENT
        chain: based on connectivity map of each atom in each residue
        metalatom:  fix1: treat as isolated atom
                    fix2: connect to donor atom (MCPB?)
        ligand: fix1: use antechamber generated prepin file to get connectivity.
                      according to https://ambermd.org/doc/prep.html the coordniate line will always start at the 11th line after 3 DUMM.
        """
        # san check
        if ligand_fix == 1 and prepi_path == None:
            raise Exception("Ligand fix 1 requires prepin_path.")
        # chain part
        for chain in self.chains_:
            for res in chain:
                for atom in res: #@shaoqz: @imp this no longer works right?
                    atom.get_connect()
        for sol in self.solvents:
            for atom in sol:
                atom.get_connect()
        # metal
        for metal in self.metalatoms:
            metal.connect = []
        if metal_fix == 1:
            pass
        if metal_fix == 2:
            raise Exception("TODO: Still working on 2 right now")

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
                            lig._find_atom_name(lp[0]).connect.append(
                                lig._find_atom_name(lp[1])
                            )
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

    def get_connectivty_table( #@shaoqz: ok seems not using
        self, ff="GAUSSIAN", metal_fix=1, ligand_fix=1, prepi_path=None
    ):
        """
        get connectivity table with atom index based on 'ff' settings:
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
        self.get_connect(metal_fix, ligand_fix, prepi_path)

        # write str in order
        # Note: Only write the connected atom with larger id
        a_id = 0
        for chain in self.chains_:
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

    def get_all_protein_atom(self): #@shaoqz: @nu
        """
        get a list of all protein atoms
        return all_P_atoms
        """
        all_P_atoms = []
        for chain in self.chains_:
            for residue in chain:
                all_P_atoms.extend(residue.atoms)
        return all_P_atoms

    def get_all_residue_unit(self, ifsolvent=0): #@shaoqz:@nu
        all_r_list = []
        for chain in self.chains_:
            for resi in chain:
                all_r_list.append(resi)

        for resi in self.ligands:
            all_r_list.append(resi)

        for resi in self.metalatoms:
            all_r_list.append(resi)

        if ifsolvent:
            for resi in self.solvents:
                all_r_list.append(resi)

        return all_r_list

    def get_residue(self, id): #@shaoqz: @imp duplicated ones?
        """
        re-search residue id with all residues count togethor from 1.
        ----------
        return a residue object
        """
        all_resi = self.get_all_residue_unit()
        for resi in all_resi:
            if resi.id == int(id):
                return resi

    def chains(self) -> List[Chain]:
        """Getter for the list of Chain() objects contained within the Structure() object."""
        self.chains_.sort(key=lambda c: c.name())
        # TODO(CJ): should this be a deep copy or no? #@shaoqz: No
        return self.chains_

    def get_atom_id(self): #@shaoqz: @nu
        """
        return a list of id of all atoms in the structure
        """
        atom_id_list = []
        for chain in self.chains_:
            for res in chain:
                for atom in res:
                    if atom.id == None:
                        raise Exception(
                            "Detected None in chain "
                            + str(chain.id)
                            + str(res.id)
                            + " "
                            + atom.name
                        )
                    atom_id_list.append(atom.id)

        for metal in self.metalatoms:
            if metal.id == None:
                raise Exception("Detected None in metal " + metal.name)
            atom_id_list.append(metal.id)

        for lig in self.ligands:
            for atom in lig:
                if atom.id == None:
                    raise Exception("Detected None in ligands", res.id, atom.name)
                atom_id_list.append(atom.id)

        for sol in self.solvents:
            for atom in sol:
                if atom.id == None:
                    raise Exception("Detected None in solvent", res.id, atom.name)
                atom_id_list.append(atom.id)

        return atom_id_list

    def get_atom_charge(self, prmtop_path): #@shaoqz: @nu
        """
        requires generate the stru using !SAME! PDB as one that generate the prmtop.
        """
        pass

    def get_atom_type(self, prmtop_path): #@shaoqz: @nu
        """
        requires generate the stru using !SAME! PDB as one that generate the prmtop.
        """
        # get type list
        with open(prmtop_path) as f:
            type_list = []
            line_index = 0

            for line in f:

                line_index = line_index + 1  # current line

                if line.strip() == r"%FLAG POINTERS":
                    format_flag = line_index
                if line.strip() == r"%FLAG AMBER_ATOM_TYPE":
                    type_flag = line_index

                if "format_flag" in dir():
                    if line_index == format_flag + 2:
                        N_atom = int(line.split()[0])
                        del format_flag

                if "type_flag" in dir():
                    if (
                        line_index >= type_flag + 2
                        and line_index <= type_flag + 1 + ceil(N_atom / 5)
                    ):
                        for i in line.strip().split():
                            type_list.append(i)
        # assign type to atom
        for chain in self.chains_:
            for res in chain:
                for atom in res:
                    atom.type = type_list[atom.id - 1]
        for atom in self.metalatoms:
            if atom.id == None:
                raise Exception("Detected None in metal " + atom.name)
            atom.type = type_list[atom.id - 1]
        for lig in self.ligands:
            for atom in lig:
                if atom.id == None:
                    raise Exception("Detected None in ligands", res.id, atom.name)
                atom.type = type_list[atom.id - 1]
        for sol in self.solvents:
            for atom in sol:
                if atom.id == None:
                    raise Exception("Detected None in solvent", res.id, atom.name)
                atom.type = type_list[atom.id - 1]

    def get_resi_dist(self, r1, r2, method="mass_center"): #@shaoqz: @nu
        """
        r1: residue 1. (residue_obj)
        r2: residue 2. (residue_obj)
        method: The method to narrow the residue down to a point.
                - mass_center
                - ...
        """
        if method == "mass_center":
            p1 = r1.get_mass_center()
            p2 = r2.get_mass_center()

        D = get_distance(p1, p2)

        return D

    def __bool__(self) -> bool:
        """Enables running assert Structure(). Checks if there is anything in the structure."""
        return bool(len(self.chains_))

    def __eq__(self, other: Structure) -> bool:
        """Comparison operator for other Structure() objects. Checks first if both have same chain names and then if each named chain is identical."""
        if set(self.chain_mapper.keys()) != set(other.chain_mapper.keys()):
            return False

        chain_name: Chain
        other_chain: Chain
        for chain_name, self_chain in self.chain_mapper.items():
            other_chain = other.chain_mapper[chain_name]
            if not self_chain.same_sequence(other_chain):
                return False
        return True #@shaoqz: so this comparsion is only in sequence level. This does not really make sense. Different levels of comparsion
                    #         is needed for a pair structure
                    #         for example:
                    #            - if the *coordinate* of every atom is the same
                    #            - if the *topology* is the same. (presence of atoms and their connectivity)
                    #            - if the *sequence* is the same. (this is basic works as the topology one as topology inside and between
                    #              residues are relatively conserved, but there are exceptions, e.g. different ligand both name LIG)

    def __ne__(self, other: Structure) -> bool:
        """Negation operator for other Structure() objects. Inverstion of Structure.__eq__()."""
        return not (self == other)

    def insert_residue(self, new_res: Residue) -> None:
        """Inserts a new Residue() object into the Structure(). If the exact Residue (chain_id, name, residue_id) already
        exists, the new Residue overwrites it. If the new Residue specifies a Chain() that does not exist, a new chain is made.
        """
        chain_name: str = new_res.chain()
        if not self.has_chain(chain_name):
            new_chain: Chain = Chain(chain_name, [new_res]) #@shaoqz: give a warning
            self.chains_.apppend(new_chain)
            self.chain_mapper[chain_name] = new_chain
        else:
            self.chain_mapper[chain_name].insert_residue(new_res)

        self.chains_ = list(self.chain_mapper.values()) #@shaoqz: why need this
        self.chains_.sort(key=lambda c: c.name()) #@shaoqz: should this be in the 1st if block?

    def get_residue(self, target_key: str) -> Union[None, Residue]: #@shaoqz: we need to that gives a reference. In the case of editing the structure.
        """Given a target_key str of the Residue() residue_key ( "chain_id.residue_name.residue_number" ) format,
        a deepcopy of the corresponding Residue() is returned, if it exists. None is returned if it cannot be found."""
        for chain in self.chains_:
            for res in chain.residues():
                if res.residue_key == target_key:
                    return deepcopy(res)
        return None

    def remove_residue(self, target_key: str) -> None:
        """Given a target_key str of the Residue() residue_key ( "chain_id.residue_name.residue_number" ) format,
        the Residue() is removed if it currently exists in one of the child Chain()'s. If the Chain() is empty after this
        removal, the chain is deleted."""
        (chain_name, _, _) = target_key.split(".") #@shaoqz: why not use this in get lolll
        if self.has_chain(chain_name):
            self.chain_mapper[chain_name].remove_residue(target_key)
            if self.chain_mapper[chain_name].empty():
                self.remove_chain(chain_name)

    def chain_names(self) -> List[str]:
        """Returns a list of all the chain names for the Structure()"""
        return list(self.chain_mapper.keys())


def compare_structures(left: Structure, right: Structure) -> Dict[str, List[str]]:
    """Compares two Structure() objects and returns a dict() of missing Residues with format:

    {'left': ['residue_key1','residue_key1',..],
     'right': ['residue_key1','residue_key1',..]
         }
    """
    result = {"left": [], "right": []}
    left_keys: Set[str] = set(left.residue_keys())
    right_keys: Set[str] = set(right.residue_keys())

    result["left"] = list(filter(lambda ll: ll not in right_keys, left_keys))
    result["right"] = list(filter(lambda rr: rr not in left_keys, right_keys))
    return result


def merge_right(left: Structure, right: Structure) -> Structure: #@shaoqz: I believe there will be a bug due to the treatment of insert with same id
    """Merges Residue() and derived objects from left Structure() to right Structure(), making sure that ALL Residue() and
        Residue() derived objects from the left are in the right. Note that the reverse is not applied and that elements initially found only in right are NOT
        merged back over to left. Also not the resulting Structure() is a deepcopy and no changes are made to the original left or right objects.

    Example:
            #TODO(CJ): add this in
    """
    struct_cpy: Structure = deepcopy(right)
    # TODO(CJ): make this a method
    left.chains_ = sorted(left.chains_, key=lambda c: c.name())
    struct_cpy.chains_ = sorted(struct_cpy.chains_, key=lambda c: c.name())
    # this is the case where there is straight up a missing chain
    for cname, chain in left.chain_mapper.items():
        if not struct_cpy.has_chain(cname):
            struct_cpy.insert_chain(deepcopy(chain))

    right_keys = struct_cpy.residue_keys() #@shaoqz: @imp2 why is this part needed as all res is stored in chain.
    for lkey in left.residue_keys():
        if lkey not in right_keys:
            struct_cpy.insert_residue(left.get_residue(lkey))
    return struct_cpy