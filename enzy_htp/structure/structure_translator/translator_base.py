"""The TranslatorBase describes the API for supporting the translation of a Structure's canonical amino acid Residue and Atom names. 
In this scheme, AmberMD naming represents the standard naming scheme. Mapping is executed in the TranslatorBase.init_mappings() function,
which will call on static, pre-defined calls to TranslatorBase.register_mapping() function. For more information and examples, see 
the rosetta version in enzy_htp/structure/translate_structure/rosetta_translator.py

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-05-27
"""

from typing import Dict, Tuple, List
from abc import ABC, abstractmethod

from plum import dispatch

from enzy_htp.core import _LOGGER

from ..atom import Atom
from ..residue import Residue
from ..structure import Structure


class TranslatorBase(ABC):
    """Defines the API for translating to/from a given naming scheme. There must be one for each
    naming scheme, the standard naming scheme is AmberMD and translation is primarily targeted at 
    canonical amino acids. All translation occurs through overloaded versions of the to_standard()
    and from_standard() methods.

    Attributes:
        TO_STANDARD: A dict() with (key, value) pairs of (translated_resname, standard_res_name) and 
            ((translated_resname, translated_atom_name), (standard_resname, standard_atom_name)) for each 
            residue code and atom name, respectively.
        FROM_STANDARD: A dict() with (key, value) pairs of (standard_res_name, translated_resname) and 
            ((standard_resname, standard_atom_name), (translated_resname, translated_atom_name)) for each
            residue code and atom name, respectively.
    """
    def __init__(self):
        """Simplistic constructor that initializes variables then calls child-instantiated init_mappings() method."""
        self.TO_STANDARD = dict()
        self.FROM_STANDARD = dict()
        self.init_mappings()

    
    def register_mapping(self,
                        s_rname:str,
                        s_atoms:List[str],
                        t_rname:str,
                        t_atoms:List[str]) -> None:
        """Workhorse method where one defines how the same residue is described in multiple 
        software packages. s_ prefix is for standard, t_ prefix is for translatred. Function works 
        for only one Residue. Note that the number of Atom()'s must be the same and that
        the supplied atom names in s_atoms and 
        
        Args:
            s_rname: What is the 3-letter Residue name in AmberMD?
            s_atoms: What are the atom names for the Residue in AmberMD?
            t_rname: What is the 3-letter Residue name in the other naming scheme?
            s_atoms: What are the atom names for the Residue in the other name scheme?

        Returns:
            Nothing.
        """

        n_s_atoms:int = len(s_atoms)
        n_t_atoms:int = len(t_atoms)
        
        if n_s_atoms != n_t_atoms:
            err_msg = f"Error! Number of standard-named and translated atoms must be the same, but got {n_s_atoms} vs {n_t_atoms}"
            _LOGGER.error(err_msg)
            raise ValueError(err_msg)

        for s_aa, t_aa in zip(s_atoms, t_atoms):
            s_key:Tuple[str, str] = (s_rname, s_aa)
            t_key:Tuple[str, str] = (t_rname, t_aa)
            self.TO_STANDARD[t_key] = s_key
            self.FROM_STANDARD[s_key] = t_key
            self.TO_STANDARD[t_rname] = s_rname
            self.FROM_STANDARD[s_rname] = t_rname
    
    @dispatch
    def to_standard(self, stru:Structure) -> None:
        for res in stru.residues:
            if not res.is_canonical():
                continue
            
            for atom in res.atoms:
                self.to_standard(atom)

            self.to_standard(res)

    @dispatch
    def to_standard(self, res:Residue) -> None:

        value = self.TO_STANDARD.get(res.name, None)

        if value is None:
            return

        res.name = value


    @dispatch
    def to_standard(self, atom:Atom) -> None:
        
        key = (atom.parent.name, atom.name)
        value = self.TO_STANDARD.get(key, None)

        if value is None:
            return

        atom.name = value[1]

    @dispatch
    def from_standard(self, stru:Structure) -> None:
        for res in stru.residues:
            if not res.is_canonical():
                continue
            
            for atom in res.atoms:
                self.from_standard(atom)

            self.from_standard(res)


    @dispatch
    def from_standard(self, res:Residue) -> None:
        
        value = self.FROM_STANDARD.get(res.name, None)

        if value is None:
            return

        res.name = value

    @dispatch
    def from_standard(self, atom:Atom) -> None:
        
        key = (atom.parent.name, atom.name)
        value = self.FROM_STANDARD.get(key, None)

        if value is None:
            return

        atom.name = value[1]


    @abstractmethod
    def naming_style(self) -> str:
        """What naming scheme is in use?"""
        pass

    @abstractmethod
    def init_mappings(self) -> None:
        """Mappings are defined here through various calls to regiset_mapping()"""
        pass

