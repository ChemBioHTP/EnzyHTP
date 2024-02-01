"""This module contains the actual ResidueCap classes that connect with free valences in the 
StructureRegion()'s. It includes the following classes:

    + ResidueCap - The abstract base class defining the API
    + HCap       - Implements a -H cap.
    + CH3Cap     - Implements a -CH3 cap.
    + COCH3Cap   - Implements a -COCH3 cap.
    + NHCH3Cap   - Implements a -NHCH3 cap.

All of these ResidueCap()'s are summarized in the below Dict[str,ResidueCap]:

    + SUPPORTED_CAPS

As well as the function that can actually apply a ResidueCap to a Residue:

    + cap_residue_free_terminal

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2024-01-31
"""
from __future__ import annotations
from abc import abstractmethod, ABC
from copy import deepcopy
import numpy as np
from scipy.spatial.transform import Rotation as R
from typing import List, Tuple, Dict

from ..structure import Structure, Residue, Atom, Chain

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode
from enzy_htp.core.math_helper import rotation_matrix_from_vectors
from enzy_htp import chemical as chem


class ResidueCap(Residue, ABC):
    """The residue cap that link to the corresponding
    Residue while not having the children atoms actually in
    that Residue. It ensures that given the defined terminal atom, it can 
    generate a realistic geometry.

    This is mainly used in handling cases like capping atoms that
    do not exist in the real structure.

    Attributes:
        link_residue: The Residue() being capped.
        link_atom: THe terminal Atom() that needs capping specifically.
        socket_atom: Atom() that connects to the link_atom in the normal Structure().
        _terminal_type: What type of terminal is it? n-term, c-term, etc.
        atoms: A List[Atom] of the Atom()'s in the ResidueCap.
    
    Properties:
        idx: The idx of the Residue() being capped.
        chain: The chain name of the Residue() being capped.
        name: The residue name of the Residue() being capped.
        term_type: The property getter of _terminal_type
        net_charge: The net charge of the ResidueCap() chemical group.
        """
    def __init__(self, link_residue: Residue,
                 link_atom: Atom, socket_atom: Atom,
                    terminal_type:str):
        """Simplistic ctor. Will throw if link_atom is NOT in link_residue."""
        self.link_residue = link_residue
        self.link_atom = link_atom
        self.socket_atom = socket_atom
        if link_atom not in link_residue.atoms:
            _LOGGER.error("link atom must be one of the link residue atoms")
            raise ValueError

        self._rtype = chem.ResidueType.RESIDUECAP
        self._terminal_type = terminal_type 

        if self.term_type == 'nterm':
            self.atoms = self.get_nterm_atoms()
        elif self.term_type == 'cterm':
            self.atoms = self.get_cterm_atoms()

        self.set_ghost_parent()

        self.anchor()


    def root(self, _count: int = 0 ):
        
        return self.link_atom.root()



    @abstractmethod
    def anchor(self) -> None:
        """Function that moves the ResidueCap()'s Atom()'s to be in the correct position for the type of Cap that it is."""
        pass

    @abstractmethod
    def get_nterm_atoms(self) -> List[Atom]:
        """Create the default Atom() objects with default names and coordinates for n-terminal verison of the Cap."""
        pass

    @abstractmethod
    def get_cterm_atoms(self) -> List[Atom]:
        """Create the default Atom() objects with default names and coordinates for c-terminal verison of the Cap."""
        pass

    @property
    @abstractmethod
    def FIXED_CHARGE_MAP(self) -> Dict[str, float]:
        """Each specialization must have a charge map.""" 
        pass

    @property
    @abstractmethod
    def CAP_BOND_DISTANCE(self) -> Dict[str, float]:
        """Each specialization must have a charge map.""" 
        pass

    @property
    @abstractmethod
    def cap_type(self) -> str:
        """What type of cap is this? Hardcoded value specific to each specialization."""
        pass


    @property
    @abstractmethod
    def net_charge(self) -> int:
        """The int() net charge of the ResidueCap()."""
        pass

    @property
    def idx(self) -> int:
        """Gets the index of the link_residue."""
        return self.link_residue.idx

    @property
    def chain(self) -> Chain:
        """Gets the chain name of the link_residue."""
        return self.link_residue.parent

    @property
    def name(self) -> str:
        """Gets the residue name of the link_residue."""
        return self.link_residue.name
    
    @property
    def term_type(self) -> str:
        """What type of capping is this? n-term, c-term, etc."""
        return self._terminal_type

    @term_type.setter
    def term_type(self, val_in:str) -> None:
        """Setter for the terminal type of the ResidueCap(). Does not perform validation."""
        self._terminal_type = val_in

    def init_fixed_atomic_charges(self) -> None:
        """init hard coded fixed atomic charges.
        Used when init_charge and fixed charge stragety is used."""
        for atom in self.atoms:
            atom.charge = self.FIXED_CHARGE_MAP[atom.name]

    def is_residue_cap(self) -> bool:
        """Does this class cap a Residue? Always True."""
        return True

    def is_nterm_cap(self) -> bool:
        """Does this ResidueCap attach to the n-terminal end of an amino acid Residue?"""
        return self.term_type == "nterm"

    def is_cterm_cap(self) -> bool:
        """Does this ResidueCap attach to the c-terminal end of an amino acid Residue?"""
        return self.term_type == "cterm"

    def __str__(self) -> str:
        """str() rep of the ResidueCap() that shows its chemical type as well as what Atom()'s it is linked to."""
        return f"{self.cap_type}Cap(link_atom={self.link_atom.name},link_residue={self.link_residue.key_str})"


class HCap(ResidueCap):
    """A ResidueCap() which uses a -H (hydrogen) as a cap. Adds a total of 1 atom1, with the 
    
    n-terminal version having atom names:       HP11 
    and c-terminal version having atom names:   HP21

    Can generally be used for both n-terminal/c-teriminal and other non-amino acid bonds
    """
    FIXED_CHARGE_MAP = {
        # C-ter
        "HP21" : 0.0,
        # N-ter
        "HP11" : 0.0,
    } # improve this
    """Fixed charge for a -H attached to the N/C-ter
    of a residue."""

    CAP_BOND_DISTANCE:Dict[str, float] = {
        "N" : 1.1,
        "C" : 1.0,
    }
    """Maps the name of the link_atom to appropriate bond distance."""

    @property
    def cap_type(self) -> str:
        """hard coded cap type"""
        return "H"

    @property
    def net_charge(self) -> str:
        """hard coded net charge"""
        return 0

    def anchor(self) -> None:
        """Straightforward anchoring that is identical for both n-terminal and c-terminal capping."""
        p1 = np.array(self.link_atom.coord)
        p2 = np.array(self.socket_atom.coord)

        direction = p2 - p1
        direction /= np.linalg.norm(direction)
    
        cap_bond_distance = self.CAP_BOND_DISTANCE[self.link_atom.name]

        self.atoms[0].coord = p1 + cap_bond_distance*direction

    def get_nterm_atoms(self) -> List[Atom]:
        """Create the default n-terminal version of the HCap with appropriate names."""
        return [Atom(name='HP11', coord=[0.0, 0.0, 0.0], element='H' )]

    def get_cterm_atoms(self) -> List[Atom]:
        """Create the default c-terminal version of the HCap with appropriate names."""
        return [Atom(name='HP21', coord=[0.0, 0.0, 0.0], element='H' )]


class CH3Cap(ResidueCap):
    """A ResidueCap() which uses a -CH3 (methyl) as a cap. Adds a total of 4 atoms, with the 
    
    n-terminal version having atom names:       CP1, HP11, HP12, HP13
    and c-terminal version having atom names:   CP2, HP21, HP22, HP23

    Can generally be used for both n-terminal/c-teriminal and other non-amino acid bonds
    """
    FIXED_CHARGE_MAP:Dict[str, float] = {
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
    """Fixed charge for a -CH3 attached to the N/C-ter
    of a residue. Modified from NME and ACE in aminon/ct12.lib"""

    CAP_BOND_DISTANCE:Dict[str, float] = {
        "N" : 1.5,
        "C" : 1.5,
    }
    """Maps the name of the link_atom to appropriate bond distance."""

    def get_nterm_atoms(self) -> List[Atom]:
        """Create the default n-terminal version of the CH3Cap with appropriate names."""
        return [
        Atom(name="CP1" , coord=[0.000,  0.000,  0.000], element='C'),
        Atom(name="HP11", coord=[0.360, -1.029,  0.000], element='H'),
        Atom(name="HP12", coord=[0.360,  0.514,  0.891], element='H'),
        Atom(name="HP13", coord=[0.360,  0.514, -0.891], element='H'),]
        
    def get_cterm_atoms(self) -> List[Atom]:
        """Create the default c-terminal version of the CH3Cap with appropriate names."""
        return [
        Atom(name="CP2" , coord=[0.000,  0.000,  0.000], element='C'),
        Atom(name="HP21", coord=[0.360, -1.029,  0.000], element='H'),
        Atom(name="HP22", coord=[0.360,  0.514,  0.891], element='H'),
        Atom(name="HP23", coord=[0.360,  0.514, -0.891], element='H'),]

    @property
    def cap_type(self) -> str:
        """What type of """
        return "CH3"

    @property
    def net_charge(self) -> str:
        """hard coded net charge"""
        return 0

    def anchor(self) -> None:
        """Straightforward anchoring that is identical for both n-terminal and c-terminal capping."""
        
        p1 = np.array(self.link_atom.coord)
        p2 = np.array(self.socket_atom.coord)

        direction = p2 - p1
        direction /= np.linalg.norm(direction)
        rot_mat = rotation_matrix_from_vectors(np.array([1,0,0]), direction)
        cap_bond_distance = self.CAP_BOND_DISTANCE[self.link_atom.name]

        for aa in self.atoms:
            pt = np.array(aa.coord)
            pt = np.transpose(np.matmul(rot_mat, np.transpose( pt  )))
            pt += (direction*cap_bond_distance) + p1 
            aa.coord = pt


class NHCH3Cap(ResidueCap):
    """A ResidueCap() which uses a -NHCH3 (methyl-amide) as a cap. Adds a total of 6 atoms, with the 
    
    n-terminal version having atom names:       NP1, HP1, CP1, HP11, HP12, HP13 
    and c-terminal version having atom names:   NP2, HP2, CP2, HP21, HP22, HP23

    Should really only be used for c-teriminal capping.
    """

    FIXED_CHARGE_MAP:Dict[str, float] = {
        # N-terminal
        'NP1': 0.0,
        'HP1': 0.0,
        'CP1': 0.0,
        'HP11': 0.0,
        'HP12': 0.0,
        'HP13': 0.0,
        # C-terminal
        "NP2": 0.0,
        "HP2": 0.0,
        "CP2": 0.0,
        "HP21": 0.0,
        "HP22": 0.0,
        "HP23": 0.0,
    }
    """Fixed charge for a NHCH3 attached to the N/C-ter of a Residue()."""


    CAP_BOND_DISTANCE:Dict[str, float] = {
        'N': 1.33,
        'C': 1.33
    }
    """Maps the name of the link_atom to appropriate bond distance."""

    def get_nterm_atoms(self) -> List[Atom]:
        """Create the default n-terminal version of the NHCH3Cap with appropriate names."""
        return [
            Atom(name='NP1',  coord=[ 0.000,   0.000,  0.000], element='N'),
            Atom(name='HP1',  coord=[ -0.5,   -0.860,  0.000], element='H'),
            Atom(name='CP1',  coord=[ 1.500,   0.000,  0.000], element='C'),
            Atom(name='HP11', coord=[ 1.864,   0.000, -1.027], element='H'),
            Atom(name='HP12', coord=[ 1.864,  -0.889,  0.515], element='H'),
            Atom(name='HP13', coord=[ 1.864,   0.889,  0.515], element='H'),
        ]


    def get_cterm_atoms(self) -> List[Atom]:
        """Create the default c-terminal version of the NHCH3Cap with appropriate names."""
        return [
            Atom(name='NP2',  coord=[ 0.000,   0.000,  0.000], element='N'),
            Atom(name='HP2',  coord=[ -0.5,   -0.860,  0.000], element='H'),
            Atom(name='CP2',  coord=[ 1.500,   0.000,  0.000], element='C'),
            Atom(name='HP21', coord=[ 1.864,   0.000, -1.027], element='H'),
            Atom(name='HP22', coord=[ 1.864,  -0.889,  0.515], element='H'),
            Atom(name='HP23', coord=[ 1.864,   0.889,  0.515], element='H'),
        ]

    def net_charge(self) -> int:
        """Net charge always zero for this neutral group."""
        return 0

    def anchor(self) -> None:
        """Straightforward anchoring that is mostly identifcal for both n-terminal and c-terminal capping. Will warn when 
        user attempts to cap with NHCH3 on n-terminal side as this creates an unusual and unrealistic NH-NH bond."""

        if self.is_nterm_cap():
            _LOGGER.warning("WARNING: You are trying to add a methylamide cap on the n-terminal side of an amino acid!")
        
        p0 = np.array(self.link_atom.coord)
        p1 = np.array(self.socket_atom.coord)
        p2 = np.array(self.socket_atom.parent.find_atom_name('CA').coord)

        direction0 = p1 - p0        
        direction0 /= np.linalg.norm(direction0)
        
        direction = p2 - p1
        direction /= np.linalg.norm(direction)
        rot_mat = rotation_matrix_from_vectors(np.array([1,0,0]), direction)
        
        cap_bond_distance = self.CAP_BOND_DISTANCE[self.link_atom.name]

        for aa in self.atoms:
            pt = np.array(aa.coord)
            pt = np.transpose(np.matmul(rot_mat, np.transpose( pt  )))
            pt += (direction0*cap_bond_distance) + p0
            aa.coord = pt

    @property
    def cap_type(self) -> str:
        """This is the NHCH3/methylamide ResidueCap."""
        return "NHCH3"


class COCH3Cap(ResidueCap):
    """A ResidueCap() which uses a -COCH3 (acetyl) as a cap. Adds a total of 6 atoms, with the 
    
    n-terminal version having atom names:       CP0, OP1, CP1, HP11, HP12, HP13 
    and c-terminal version having atom names:   CP2, OP2, CP2, HP21, HP22, HP23

    Should really only be used for n-teriminal capping.
    """
    FIXED_CHARGE_MAP:Dict[str,float] = {
        # n-terminal
        'CP0': 0.0,
        'OP1': 0.0,
        'CP1': 0.0,
        'HP11': 0.0,
        'HP12': 0.0,
        'HP13': 0.0,
        # c-terminal
        'CP2': 0.0,
        'OP2': 0.0,
        'CP3': 0.0,
        'HP31': 0.0,
        'HP32': 0.0,
        'HP33': 0.0, 
    }
    """Fixed charge for a NHCH3 attached to the N/C-ter of a Residue()."""

    CAP_BOND_DISTANCE:Dict[str, float] = {
        "N" : 1.33,
        "C" : 1.33,
    }
    """Maps the name of the link_atom to appropriate bond distance."""

    @property
    def cap_type(self) -> str:
        """This is the COCH3/acetyl ResidueCap."""
        return "COCH3"


    def anchor(self) -> None:
        """Straightforward anchoring that is mostly identifcal for both n-terminal and c-terminal capping. Will warn when 
        user attempts to cap with COCH3 on c-terminal side as this creates an unusual and unrealistic CO-CO bond."""

        if self.is_cterm_cap():
            _LOGGER.warning("WARNING: You are trying to add an acetyl cap on the c-terminal side of an amino acid!")
            self.atoms[0].coord = self.socket_atom.parent.find_atom_name('N').coord
            self.atoms[1].coord =  self.socket_atom.parent.find_atom_name('H').coord
            direction =  np.array(self.atoms[1].coord) - np.array(self.atoms[0].coord)
            direction /= np.linalg.norm(direction)
            self.atoms[1].coord = 1.2*direction + np.array(self.atoms[0].coord)
            p0 = np.array(self.link_atom.coord)
            p1 = np.array(self.socket_atom.coord)
            p2 = np.array(self.socket_atom.parent.find_atom_name('CA').coord)

            direction0 = p1 - p0        
            direction0 /= np.linalg.norm(direction0)
            
            direction = p2 - p1
            direction /= np.linalg.norm(direction)
            rot_mat = rotation_matrix_from_vectors(np.array([1,0,0]), direction)
            
            cap_bond_distance = self.CAP_BOND_DISTANCE[self.link_atom.name]

            for aa in self.atoms[2:]:
                pt = np.array(aa.coord)
                pt = np.transpose(np.matmul(rot_mat, np.transpose( pt  )))
                pt += (direction0*cap_bond_distance) + p0
                aa.coord = pt

        elif self.is_nterm_cap():
            self.atoms[0].coord = self.socket_atom.parent.find_atom_name('C').coord
            self.atoms[1].coord = self.socket_atom.parent.find_atom_name('O').coord
   
            p0 = np.array(self.socket_atom.parent.find_atom_name('C').coord)
            p1 = np.array(self.socket_atom.parent.find_atom_name('CA').coord)

            direction0 = p1 - p0
            direction0 /= np.linalg.norm(direction0)
            rot_mat = rotation_matrix_from_vectors(np.array([1,0,0]), direction0)
            cap_bond_distance = self.CAP_BOND_DISTANCE[self.link_atom.name]

            for aa in self.atoms[2:]:
                pt = np.array(aa.coord)
                pt = np.transpose(np.matmul(rot_mat, np.transpose( pt  )))
                pt += p0
                aa.coord = pt


    def net_charge(self) -> int:
        """Net charge always zero for this neutral group."""
        return 0

    def get_nterm_atoms(self) -> List[Atom]:
        """Create the default n-terminal version of the COCH3Cap with appropriate names."""
        return [
            Atom(name='CP0',  coord=[ 0.000 ,   0.000,   0.000], element= 'C'),
            Atom(name='OP1',  coord=[-0.610 ,   1.056,   0.000], element= 'O'),
            Atom(name='CP1',  coord=[ 1.507 ,   0.000,   0.000], element= 'C'),
            Atom(name='HP11', coord=[ 1.871 ,   0.000,   1.028], element= 'H'),
            Atom(name='HP12', coord=[ 1.871 ,  -0.890,  -0.514], element= 'H'),
            Atom(name='HP13', coord=[ 1.871 ,   0.890,  -0.514], element= 'H'),
        ]
    
    def get_cterm_atoms(self) -> List[Atom]:
        """Create the default c-terminal version of the COCH3Cap with appropriate names."""
        return [
            Atom(name='CP2',  coord=[ 0.000 ,   0.000,   0.000], element= 'C'),
            Atom(name='OP2',  coord=[-0.610 ,   1.056,   0.000], element= 'O'),
            Atom(name='CP3',  coord=[ 1.507 ,   0.000,   0.000], element= 'C'),
            Atom(name='HP31', coord=[ 1.871 ,   0.000,   1.028], element= 'H'),
            Atom(name='HP32', coord=[ 1.871 ,  -0.890,  -0.514], element= 'H'),
            Atom(name='HP33', coord=[ 1.871 ,   0.890,  -0.514], element= 'H'),
        ]


def cap_residue_free_terminal(
        res: Residue, cap_name: str, terminal_type: str) -> List[Atom]:
    """Cap the free terminal of {terminal_type} of the give residue {res}
    using cap of {cap_name}. Implementatino function NOT TO BE CALLED BY USERS DIRECTLY. 
    Args:
        res:
            the target residue
        cap_name:
            the name of the cap
        terminal_type:
            the free terminal type.
            options: [nterm, cterm]
    Returns:
        The List[Atom] of capping atoms.
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
    
    residue_cap:ResidueCap = SUPPORTED_CAPS[cap_name](res, start_atom, end_atom, terminal_type)
    
    return residue_cap.atoms

SUPPORTED_CAPS: Dict[str, callable] = {
    # chemical formulas first
    "CH3" : CH3Cap,
    "H" : HCap,
    "NHCH3": NHCH3Cap,
    "COCH3": COCH3Cap,
    # common names for redundancy
    "methyl": CH3Cap,
    "hydrogen": HCap,
    "methylamide": NHCH3Cap,
    "acetate": COCH3Cap
}
"""Mapper that converts a cap name into a constructor. Both the chemical formula and common names for each method are provided."""
