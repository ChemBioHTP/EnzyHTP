"""TODO DOCUMENATION"""
# TODO update key when chain's name is changed
from __future__ import annotations

import numpy as np
from typing import Tuple, List
from plum import dispatch
from .atom import Atom
from enzy_htp.chemical import ResidueType, METAL_MAPPER, THREE_LETTER_AA_MAPPER, RD_SOLVENT_LIST, RD_NON_LIGAND_LIST


class Residue:
    """TODO DOCUMENTATION"""
    def __init__(self, residue_key, atoms):
        self.atoms = atoms
        self.residue_key = residue_key
        (chain, name, num) = self.residue_key.split('.')
        self.chain_ = chain
        self.name = name
        self.num_ = int(num)
        self.rtype_ = ResidueType.UNKNOWN
        line_idxs = np.array(list(map(lambda aa: aa.line_idx, self.atoms)))
        self.min_line_ = np.min(line_idxs)
        self.max_line_ = np.max(line_idxs)
        # TODO what are the checks that we should be doing here?

    def atom_list(self) -> List[Atom]:
        return self.atoms

    def __determine_residue_type( self ):
        # TODO finish this algorithm
        # 1. canoncial code => canonical
        # 2. Solvent = WAT or HOH, also solvent ions (NA+, CL-)
        # 3. Metal-center: only from metal center residue 
        #.4. Non-Canonical/Ligand => similar but ligand will be in its own chain... ligand will NEVER be in the same cahin as a canonical amino acid
        pass

    def empty_chain(self) -> bool:
        return not len(self.chain_.strip())

    def chain(self) -> str:
        return self.chain_

    def set_chain(self, val : str ) -> None:
        # TODO:CJ maybe fix this? its kinda janky
        for idx in range(len(self.atoms)):
            self.atoms[idx].chain_id = val
        self.chain_ = val

    def num(self) -> int:
        return self.num_

    def min_line(self) -> int:
        return self.min_line_
    
    def max_line(self) -> int:
        return self.max_line_

    def is_metal(self) -> bool:
        return self.name in METAL_MAPPER

    def is_canonical(self) -> bool:
        return self.name in THREE_LETTER_AA_MAPPER

    def line_range(self) -> Tuple[int,int]:
        return (self.min_line(), self.max_line())

    def neighbors(self, other : Residue) -> bool:
        """Checks if two residues from the same PDB are neighbors."""
        return (abs(self.min_line()-other.max_line()) == 1 ) or ( abs(self.max_line() - other.min_line()) == 1 )

    def is_rd_solvent(self) -> bool:
        return self.name in RD_SOLVENT_LIST

    def is_rd_non_ligand(self) -> bool:
        return self.name in RD_NON_LIGAND_LIST

    @dispatch
    def rtype(self, new_rtype: ResidueType) -> None:
        pass

    @dispatch
    def rtype( self ) -> ResidueType:
        pass
    #TODO add operator overloading to move the chain?


#class Residue(Child):
#    """
#    -------------
#    initilize from
#    PDB:        Residue.fromPDB(resi_input, resi_id=pdb_l.resi_id, input_type='PDB_line' or 'line_str' or 'file' or 'path')
#    raw data:   Residue(atoms, resi_id, resi_name)
#    -------------
#    id
#    name
#    parent # the whole chain
#
#    atoms = [atom_obj, ...]
#
#    # --after metalatom.get_donor_resi--
#    d_atom (donor atom when work as a ligand)
#    a_metal (acceptor metal)
#    
#    #TODO
#    if_art_resi
#    -------------
#    Method
#    -------------
#    set_parent
#    if_art_resi
#    deprotonate
#    _find_atom_name
#    -------------
#    __getitem__
#        Residue_obj[int]: Residue_obj.residues[int]    
#    __getattr__
#        Residue_obj.123 = Residue_obj.atoms[123-1] // index mimic (start from 1)
#        Residue_obj.CA = Residue_obj.find_atom_name('CA') // search mimic
#    __delitem__
#        del obj[int] --> obj.child[int].remove() // delete by index (start from 0)
#        del obj[str] --> obj._del_child_name(str) // delete by value
#        del obj[child] --> obj.child_list.remove(child) // delete by value
#    __len__
#        len(obj) = len(obj.child_list)
#    """
#
#    """
#    ====
#    init
#    ====
#    """
#
#    def __init__(self, atoms, resi_id, resi_name, parent=None):
#        """
#        Common part of init methods: direct from data objects
#        No parent by default. Add parent by action
#        """
#        # set parent to None
#        Child.__init__(self)
#        # add parent if provided
#        if parent != None:
#            self.set_parent(parent)
#        # adapt some children
#        self.atoms = []
#        for i in atoms:
#            i.set_parent(self)
#            self.atoms.append(i)
#        # set id
#        self.id = resi_id
#        self.name = resi_name
#
#        # clean
#        self.d_atom = None
#
#    @classmethod
#    def fromPDB(cls, resi_input, resi_id=None, input_type="PDB_line"):
#        """
#        generate resi from PDB. Require 'ATOM' and 'HETATM' lines.
#        ---------
#        resi_input = PDB_line (or line_str or file or path)
#        resi_id : int (use the number in the line by default // support customize)
#        Use PDB_line in the list to init each atom
#        """
#
#        # adapt general input // converge to a list of PDB_line (resi_lines)
#        if input_type == "path":
#            f = open(resi_input)
#            resi_lines = PDB_line.fromlines(f.read())
#            f.close()
#        if input_type == "file":
#            resi_lines = PDB_line.fromlines(resi_input.read())
#        if input_type == "line_str":
#            resi_lines = PDB_line.fromlines(resi_input)
#        if input_type == "PDB_line":
#            resi_lines = resi_input
#
#        # clean lines
#        for i in range(len(resi_lines) - 1, -1, -1):
#            if (
#                resi_lines[i].line_type != "ATOM"
#                and resi_lines[i].line_type != "HETATM"
#            ):
#                if Config.debug > 1:
#                    print("Residue.fromPDB: delete error line in input.")
#                del resi_lines[i]
#
#        # Default resi_id
#        if resi_id is None:
#            resi_id = resi_lines[0].resi_id
#        # get name from first line
#        resi_name = resi_lines[0].resi_name
#        # get child atoms
#        atoms = []
#        for pdb_l in resi_lines:
#            atoms.append(Atom.fromPDB(pdb_l))
#
#        return cls(atoms, resi_id, resi_name)
#
#    """
#    ====
#    Method
#    ====
#    """
#
#    def if_art_resi(self):
#        pass
#
#    def deprotonate(self, T_atom=None, HIP="HIE"):
#        """
#        check current protonation state.
#        deprotonate if applicable (HIP -> HIE by default)
#        ---------
#        base on T_atom if provided. (refine in the future)
#        """
#        if T_atom == None:
#            if self.name != "HIP":
#                depro_info = DeProton_map[self.name]
#            else:
#                if HIP == "HIE":
#                    depro_info = DeProton_map[self.name][0]
#                if HIP == "HID":
#                    depro_info = DeProton_map[self.name][1]
#
#            depro_resi = depro_info[0]
#            depro_atom = depro_info[1]
#
#            # delete the proton
#            del self[depro_atom]
#
#            # change the name
#            self.name = depro_resi
#
#        else:
#            # assign target atom
#            # only affect operations on residues with differences between donor atom and potential deprotonation atom (HIP HIE HID and ARG)
#            if self.name == "HIP":
#                if T_atom.name == "ND1":
#                    del self["HD1"]
#                    self.name = "HIE"
#                    return
#                if T_atom.name == "NE2":
#                    del self["HE2"]
#                    self.name = "HID"
#                    return
#
#            if self.name == "HIE":
#                if T_atom.name == "ND1":
#                    return
#                if T_atom.name == "NE2":
#                    del self["HE2"]
#                    self.name = "HID"
#                    # self.add_H('ND1')
#                    # let leap auto complete by now
#                    return
#
#            if self.name == "HID":
#                if T_atom.name == "ND1":
#                    del self["HD1"]
#                    self.name = "HIE"
#                    # self.add_H('NE2')
#                    # let leap auto complete by now
#                    return
#                if T_atom.name == "NE2":
#                    return
#
#            # find the right H to delete for ARG and LYS
#            if self.name == "ARG":
#                raise Exception(
#                    "ARG detected as donor: " + self.chain.id + str(self.id)
#                )
#
#            if self.name == "LYS":
#                # the deleted H has to be HZ1 in name
#                D1 = np.linalg.norm(
#                    np.array(self.HZ1.coord) - np.array(self.a_metal.coord)
#                )
#                D2 = np.linalg.norm(
#                    np.array(self.HZ2.coord) - np.array(self.a_metal.coord)
#                )
#                D3 = np.linalg.norm(
#                    np.array(self.HZ3.coord) - np.array(self.a_metal.coord)
#                )
#                x = min(D1, D2, D3)
#                if x == D1:
#                    del self["HZ1"]
#                    self.name = "LYN"
#                    return
#                if x == D2:
#                    del self["HZ2"]
#                    self.HZ1.name = "HZ2"
#                    self.name = "LYN"
#                    return
#                if x == D3:
#                    del self["HZ3"]
#                    self.HZ1.name = "HZ3"
#                    self.name = "LYN"
#                    return
#                raise Exception
#
#            depro_info = DeProton_map[self.name]
#            depro_resi = depro_info[0]
#            depro_atom = depro_info[1]
#            # delete the proton
#            del self[depro_atom]
#            # change the name
#            self.name = depro_resi
#
#    def add_H(self, T_atom):
#        """
#        add H on the corresponding T_atom.
#        1. make the H
#            find H name
#            find H coordinate
#        2. add H to the residue
#            use self.add()
#        """
#        pass
#
#    def rot_proton(self, T_atom):
#        """
#        rotate the dihedral relate to the target H-atom bond
#        -----
#        TODO  
#        """
#
#        if self.name == "TRP":
#            raise Exception("Error: TRP detected as donor!!!")
#
#        protons = T_atom.get_protons()
#        lp_infos = T_atom.get_lp_infos()
#        # rotate to lp direction if proton on T_atom
#        if len(protons) != 0:
#            for proton in protons:
#                bond_end1 = T_atom.get_bond_end_atom()
#                bond_end2 = bond_end1.get_bond_end_atom()
#                proton.set_byDihedral(
#                    T_atom, bond_end1, bond_end2, value=lp_infos[0]["D"]
#                )
#
#    def ifDeProton(self):
#        """
#        check if this residue add or minus proton in a pH range of 1-14. (ambiguous protonation state, potential deprotonation)
#        -------
#        base on self.name. return a bool
#        """
#        return self.name in DeProton_map.keys()
#
#    def _find_atom_name(self, name: str):
#        """
#        find atom according to the name (should find only one atom)
#        return the atom (! assume the uniqueness of name)
#        """
#        out_list = []
#        for atom in self.atoms:
#            if atom.name == name:
#                out_list.append(atom)
#        if len(out_list) > 1:
#            print(
#                "\033[32;0mShould there be same atom name in residue +"
#                + self.name
#                + str(self.id)
#                + "?+\033[0m"
#            )
#            raise Exception
#        else:
#            return out_list[0]
#
#    def _del_atom_name(self, name: str):
#        """
#        find atoms according to the name
#        delete found atoms
#        """
#        for i in range(len(self.atoms) - 1, -1, -1):
#            if self.atoms[i].name == name:
#                del self.atoms[i]
#
#    def add(self, obj):
#        """
#        1. judge obj type
#        2. set parent, add to corresponding list
#        --------
#        Do not take any id: add to the end of the residue atom list by default
#        Do not sort: the order of atom is pointless by far.
#        """
#        # list
#        if type(obj) == list:
#
#            obj_ele = obj[0]
#
#            if type(obj_ele) != Atom:
#                raise TypeError("residue.Add() method only take Atom")
#
#            obj.set_parent(self)
#            self.atoms.extend(obj)
#
#        # single building block
#        else:
#            if type(obj) != Atom:
#                raise TypeError("residue.Add() method only take Atom")
#
#            obj.set_parent(self)
#            self.atoms.append(obj)
#
#    def sort(self):
#        """
#        assign atom id based on current order in the list
#        * start from 1 in a residue
#        """
#        for index, atom in enumerate(self.atoms):
#            atom.id = index + 1
#
#    def get_mass_center(self):
#        """
#        get mass center of current residue
#        """
#        M_t = 0
#        M_cw = np.array((0.0, 0.0, 0.0))
#        for atom in self:
#            atom.get_ele()
#            a_mass = Ele_mass_map[atom.ele]
#            M_t += a_mass
#            M_cw += a_mass * np.array(atom.coord)
#
#        M_center = tuple(M_cw / M_t)
#
#        return M_center
#
#    """
#    ====
#    Special Method
#    ====
#    """
#
#    def __getitem__(self, key: int):
#        """
#        Residue_obj[int]: Residue_obj.residues[int]
#        -----
#        use residue index within the chain, start from 0
#        """
#        return self.atoms[key]
#
#    def __getattr__(self, key):
#        """
#        Residue_obj.i123 = Residue_obj.atoms[123-1] // index mimic (start from 1)
#        Residue_obj.CA = Residue_obj.find_atom_name('CA') // search mimic
#        """
#        if key == "chain":
#            return self.parent
#        # judge if a digit str, since a str will always be passed
#        if key[0] == "i":
#            key = int(key[1:])
#            return self.atoms[key - 1]
#        else:
#            return self._find_atom_name(key)
#
#    def __delitem__(self, key):
#        """
#        del obj[int] --> obj.child[int].remove() // delete by index (start from 0)
#        del obj[str] --> obj._del_child_name(str) // delete by value
#        del obj[child] --> obj.child_list.remove(child) // delete by value
#        """
#        if type(key) == int:
#            del self.atoms[key]
#        if type(key) == str:
#            return self._del_atom_name(key)
#        if type(key) == Atom:
#            self.atoms.remove(key)
#
#    def __len__(self):
#        """
#        len(obj) = len(obj.child_list)
#        """
#        return len(self.atoms)
#
#

