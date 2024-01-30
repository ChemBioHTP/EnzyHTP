"""This module handles free valence capping for StructureRegion 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2024-1-1"""
from __future__ import annotations
from copy import deepcopy
import numpy as np
from scipy.spatial.transform import Rotation as R
from typing import List, Tuple, Dict

from ..structure import Structure, Residue, Atom

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode
from enzy_htp.core.math_helper import rotation_matrix_from_vectors
from enzy_htp import chemical as chem

#TODO(CJ): add methylamide -> cterm, acetate -> nterm

# region == APIs ==

def capping_with_residue_terminals(raw_region,
                                   nterm_cap:str=None,
                                   cterm_cap:str=None,
                                   return_copy: bool= False,
                                   **kwargs):
    """cap the raw region composed by whole residues
    only. In this case, capping only on the terminal
    of the residue. (keyword: res_ter_cap)
    Args:
        raw_region: the target StructureRegion object
        nterm_cap: the cap added to the N-ter
        cterm_cap: the cap added to the C-ter
        return_copy: whether edit in place or return an capped copy
    Return:
        (if not return_copu, change the region in place)
        capped_region: the StructureRegion after capping."""
    # init
    if nterm_cap is None:
        nterm_cap = "CH3"

    if cterm_cap is None:
        cterm_cap = "CH3"

    # san check (only allow whole residues in the region)
    if not raw_region.is_whole_residue_only():
        _LOGGER.error(f"Found incomplete residue in region. ({raw_region}) "
                      "residue terminal capping only support capping for a region"
                      " composed by whole residues only.")
        raise ValueError

    free_res_mapper = raw_region.involved_residues_with_free_terminal()
    free_nter_res = free_res_mapper["n_ter"]
    free_cter_res = free_res_mapper["c_ter"]
    
    add_atoms = []
    for res in free_nter_res:
        cap_atoms = cap_residue_free_terminal(res, nterm_cap, "nterm")
        add_atoms.extend(cap_atoms)

    for res in free_cter_res:
        cap_atoms = cap_residue_free_terminal(res, cterm_cap, "cterm")
        add_atoms.extend(cap_atoms)
    
    if return_copy:
        new_region = raw_region.clone()
        new_region.atoms.extend(add_atoms)
        return new_region
    else:
        raw_region.atoms.extend(add_atoms)

# endregion

class ResidueCap(Residue):
    """the residue cap that link to the corresponding
    Residue while not having the children atoms actually in
    that Residue.
    This is mainly used in handling cases like capping atoms that
    does not exists in the real structure.
    Attribute:
        link_residue
        children/atoms
        parent = None
        link_atom:
            the linking atom in original res
        socket_atom:
            the socket from original stru to pluging the cap.
            the cap is replacing the atom.
        plug_atom:
            the pluging in atom in the cap
    
    Property:
        cap_type"""
    def __init__(self, link_residue: Residue,
                 link_atom: Atom, socket_atom: Atom,
                 atoms: List[Atom], plug_atom: Atom):
        self.link_residue = link_residue
        self.link_atom = link_atom
        self.socket_atom = socket_atom
        if link_atom not in link_residue.atoms:
            _LOGGER.error("link atom must be one of the link residue atoms")
            raise ValueError
        self.atoms = atoms
        self.plug_atom = plug_atom
        if plug_atom not in self.atoms:
            _LOGGER.error("plug atom must be one of the atoms")
            raise ValueError

        self._idx = None
        self._name = None
        self._rtype = chem.ResidueType.RESIDUECAP
        self.set_ghost_parent()
        self.set_children(atoms)
    
    @property
    def cap_type(self) -> str:
        """hard coded cap type"""
        _LOGGER.error("using abstract method.")
        raise AttributeError

    @property
    def cap_bond_vector(self) -> np.array:
        """return a vector that defines the direction of the capping bond"""
        result = np.array(self.socket_atom.coord) - np.array(self.link_atom.coord)
        result /= np.linalg.norm(result)
        return result

    def init_fixed_atomic_charges(self) -> str:
        """init hard coded fixed atomic charges.
        Used when init_charge and fixed charge stragety is used."""
        for atom in self.atoms:
            atom.charge = self.FIXED_CHARGE_MAP[atom.name]

    def is_residue_cap(self) -> bool:
        """hard coded yes"""
        return True
    
    def restore_default_coord(self):
        for atom in self.atoms:
            atom.coord = self.DEFAULT_COORDS[atom.name]

    def apply_to_geom(self, geom: Structure) -> ResidueCap:
        """make a copy of the cap for geom."""
        geom_cap = deepcopy(self)
        geom_cap.link_residue = geom.residue_mapper[self.link_residue.key()]    
        geom_cap.link_atom = geom_cap.link_residue.find_atom_name(self.link_atom.name)
        geom_cap.socket_atom = geom.get(self.socket_atom.key)
        geom_cap.restore_default_coord()
        # align to the linked residue TODO unify this. also unify the 1,0,0 state somewhere
        direction = geom_cap.cap_bond_vector
        p1 = np.array(geom_cap.link_atom.coord)
        rot_mat = rotation_matrix_from_vectors(np.array([1,0,0]), direction) # TODO is this only works for 1,0,0??
        bond_distance = self.CAP_BOND_DISTANCE[geom_cap.link_atom.name]
        for aa in geom_cap.atoms:
            pt = np.array(aa.coord)
            # rotation
            pt = np.transpose(np.matmul(rot_mat, np.transpose( pt  )))
            aa.coord = pt
            # transform
            pt += (direction*bond_distance) + p1

        return geom_cap


class CH3Cap(ResidueCap):
    """the CH3 cap
    TODO move the default atom coord in?
    TODO move bond distance in."""
    FIXED_CHARGE_MAP = {
        # C-ter
        "CP2" : -0.292800,
        "HP21" : 0.097600,
        "HP22" : 0.097600,
        "HP23" : 0.097600,
        # N-ter
        "CP1" : -0.336900,
        "HP11" : 0.112300,
        "HP12" : 0.112300,
        "HP13" : 0.112300,
    } # TODO improve this for non ter.
    """fixed charge for a -CH3 attached to the N/C-ter
    of a residue. Modified from NME and ACE in aminon/ct12.lib"""

    CAP_BOND_DISTANCE = {
        "N" : 1.5,
        "C" : 1.5,
    }
    """map the name of the link_atom to bond distance"""

    DEFAULT_COORDS = {
        "CP2"  : (0.000,  0.000,  0.000),
        "HP21" : (0.360, -1.029,  0.000),
        "HP22" : (0.360,  0.514,  0.891),
        "HP23" : (0.360,  0.514, -0.891),
        "CP1"  : (0.000,  0.000,  0.000),
        "HP11" : (0.360, -1.029,  0.000),
        "HP12" : (0.360,  0.514,  0.891),
        "HP13" : (0.360,  0.514, -0.891),
    }
    """TODO improve this"""

    @property
    def cap_type(self) -> str:
        """hard coded cap type"""
        return "CH3"

    @property
    def net_charge(self) -> str:
        """hard coded net charge"""
        return 0


class HCap(ResidueCap):
    """the H cap
    TODO move the default atom coord in?
    TODO move bond distance in."""
    FIXED_CHARGE_MAP = {
        # C-ter
        "HP21" : 0.0,
        # N-ter
        "HP11" : 0.0,
    } # improve this
    """fixed charge for a -H attached to the N/C-ter
    of a residue."""

    CAP_BOND_DISTANCE = {
        "N" : 1.1,
        "C" : 1.0,
    }
    """map the name of the link_atom to bond distance"""

    @property
    def cap_type(self) -> str:
        """hard coded cap type"""
        return "CH3"

    @property
    def net_charge(self) -> str:
        """hard coded net charge"""
        return 0


def cap_residue_free_terminal(
        res: Residue, cap_name: str, terminal_type: str) -> List[Atom]:
    """cap the free terminal of {terminal_type} of the give residue {res}
    using cap of {cap_name}
    Args:
        res:
            the target residue
        cap_name:
            the name of the cap
        terminal_type:
            the free terminal type.
            options: [nterm, cterm]
    Returns:
        a list of capping atoms
    """
    if terminal_type == 'nterm':
        start_atom:Atom = res.find_atom_name("N")
        end_atom:Atom = res.n_side_residue().find_atom_name("C")
    elif terminal_type == 'cterm':
        start_atom:Atom = res.find_atom_name("C")
        end_atom:Atom = res.c_side_residue().find_atom_name("N")
    else:
        _LOGGER.error(f"terminal_type can only be nterm or cterm. Given {terminal_type}")
        raise ValueError

    if cap_name not in SUPPORTED_CAPS:
        _LOGGER.error(f"cap_name {cap_name} not supported. Supported name: {SUPPORTED_CAPS.keys()}")
        raise ValueError
    cap_atoms = SUPPORTED_CAPS[cap_name](
        terminal_type, res, start_atom, end_atom).atoms
    cap_bond_distance = get_bond_distance(cap_name, terminal_type)

    # TODO put these in core and strutcure and in cap classes
    # TODO make this a method in cap classes
    p1 = np.array(start_atom.coord)
    p2 = np.array(end_atom.coord)
    direction = p2 - p1
    direction /= np.linalg.norm(direction)
    rot_mat = rotation_matrix_from_vectors(np.array([1,0,0]), direction)
    for aa in cap_atoms:
        pt = np.array(aa.coord)
        pt = np.transpose(np.matmul(rot_mat, np.transpose( pt  )))
        pt += (direction*cap_bond_distance) + p1 
        aa.coord = pt

    return cap_atoms

def get_bond_distance(strategy:str, side:str) -> float:
    """ """
    #TODO(CJ): update this: need methylamide and acetyl
    _BOND_DISTANCE_MAPPER:Dict[Tuple[str,str], float] = {
        ('H', 'nterm'):1.1,
        ('H', 'cterm'):1.0,
        ('CH3', 'nterm'):1.5,
        ('CH3', 'cterm'):1.5,
    }


    result:float = _BOND_DISTANCE_MAPPER.get((strategy, side), None)

    if result is None:
        _LOGGER.error(f"No bond distance for {strategy}-capping strategy on {side} side. Supported methods include:")
        for (strat, side_name) in _BOND_DISTANCE_MAPPER.keys():
            _LOGGER.error(f"strategy: {strat}, side: {side_name}")
        _LOGGER.error("Exiting...")
        exit( 1 )

    return result

def get_h_cap(end: str, res: Residue, 
              link_atom: Atom, socket_atom: Atom) -> HCap:
    """Args:
    end:
        the pattern of the end it conencts
        This determine the atom name
    res:
        the original residue.
    link_atom:
        the linking atom in original res
    socket_atom:
        the socket from original stru to pluging the cap
    
    Returns:
        a list of Atom()s in the cap"""
    if end == 'nterm':
        aname = 'HP11'
    elif end == 'cterm':
        aname = 'HP21'
    atoms = [
        Atom.from_biopandas(
            {'atom_name': aname,
            'x_coord':0.000, 'y_coord': 0.000, 'z_coord':  0.000,
            'element_symbol': "H"},
        )
    ]

    return HCap(
        link_residue=res,
        link_atom=link_atom,
        socket_atom = socket_atom,
        plug_atom = atoms[0],
        atoms = atoms,
    )

def get_ch3_cap(end: str, res: Residue,
                link_atom: Atom, socket_atom: Atom) -> CH3Cap:
    """Args:
    end:
        the pattern of the end it conencts
        This determine the atom name
    res:
        the original residue.
    link_atom:
        the linking atom in original res
    socket_atom:
        the socket from original stru to pluging the cap
    
    Returns:
        a list of Atom()s in the cap"""
    if end == 'nterm':
        anames = "CP1 HP11 HP12 HP13".split()
    elif end == 'cterm':
        anames = "CP2 HP21 HP22 HP23".split()
    atoms = [
        Atom.from_biopandas({'atom_name': anames[0], 'x_coord':0.000, 'y_coord': 0.000, 'z_coord':  0.000,'element_symbol': "C"}),
        Atom.from_biopandas({'atom_name': anames[1], 'x_coord':0.360, 'y_coord':-1.029, 'z_coord':  0.000,'element_symbol': "H"}),
        Atom.from_biopandas({'atom_name': anames[2], 'x_coord':0.360, 'y_coord': 0.514, 'z_coord':  0.891,'element_symbol': "H"}),
        Atom.from_biopandas({'atom_name': anames[3], 'x_coord':0.360, 'y_coord': 0.514, 'z_coord': -0.891,'element_symbol': "H"})
    ]
    return CH3Cap(
        link_residue=res,
        link_atom=link_atom,
        socket_atom = socket_atom,
        plug_atom = atoms[0],
        atoms = atoms,
    )

SUPPORTED_CAPS: Dict[str, callable] = {
    "CH3" : get_ch3_cap,
    "H" : get_h_cap,
}
