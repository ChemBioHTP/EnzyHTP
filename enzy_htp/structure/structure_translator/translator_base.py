from typing import Dict, Tuple, List
from abc import ABC, abstractmethod


from plum import dispatch

from enzy_htp.core import _LOGGER

from ..atom import Atom
from ..residue import Residue
from ..structure import Structure


class TranslatorBase(ABC):

    def __init__(self):
        self.TO_STANDARD = dict()
        self.FROM_STANDARD = dict()
        self.init_mappings()

    
    def register_mapping(self,
                        s_rname:str,
                        s_atoms:List[str],
                        t_rname:str,
                        t_atoms:List[str]):

        assert len(s_atoms) == len(t_atoms)

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
        pass

    @abstractmethod
    def init_mappings(self) -> None:
        pass

