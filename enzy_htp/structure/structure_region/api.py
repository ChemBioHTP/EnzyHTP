#TODO(CJ):


import numpy as np


from typing import List, Tuple, Union


from ..structure import Structure, Residue, Atom

from .capping import (
    needs_nterm_capping,
    needs_cterm_capping,
    cap_residue
    )


class StructureRegion:

    def __init__(self, atoms: List[Atom]):
        
        self.atoms_ = atoms
        self.atom_mapper_ = dict()
        for aa in self.atoms_: 
            self.atom_mapper_[aa.key] = aa

    @property
    def atoms(self) -> List[Atom]:
        return self.atoms_


    def has_atoms(self, atoms:List[Atom] ) -> bool:
        for aa in atoms:
            if not self.has_atom( aa ):
                return False
        return True
        
    def has_atom(self, atom:Atom) -> bool:
        return atom.key in self.atom_mapper_

    def get_atom(self, key:str ) -> Atom:
        #NOTE(CJ): need this part for the 
        return self.atom_mapper_[key]

    def backbone_indices(self, index:int=0):
        indices:List[str] = list()
        for aidx,aa in enumerate(self.atoms):
            if not aa.parent.is_canonical():
                continue
            
            if aa.name in "N CA C CP1 CP2".split(): #TODO(CJ): 
                indices.append( aidx + index )

        return indices

    def get_atom_index(self, atom:Atom, indexing:int=0) -> int:
        target_key:str = atom.key
        for aidx, aa in enumerate(self.atoms):
            if aa.key == target_key:
                return aidx+indexing
        else:
            assert False #TODO(CJ): better error logging here

    def closest_heavy_atom(self, atom:Atom) -> Atom:
        distances = list()
        for aa in self.atoms:
            if aa.element == 'H' or aa.parent != atom.parent:
                distances.append(1000)
            else:
                distances.append( atom.distance_to(aa) )
        return self.atoms[np.argmin(distances)]


    def __getitem__(self, key:int) -> Atom:
        return self.atoms_[key]

def create_structure_region(stru:Structure,
                            residue_list:List[Tuple[str,int]],
                            nterm_cap:str=None,
                            cterm_cap:str=None
                            ) -> StructureRegion:
    """TODO(CJ)"""

    if nterm_cap is None:
        nterm_cap = "CH3"

    if cterm_cap is None:
        cterm_cap = "CH3"

    atoms:List[Atom] = list()

    for rl in residue_list:

        res:Residue=stru.get_residue(f"{rl[0]}.{rl[1]}")

        residue_atoms:List[Atom] = list()

        if needs_nterm_capping(res, stru, residue_list):
            residue_atoms.extend(cap_residue(res, stru, nterm_cap, 'nterm'))

        residue_atoms.extend( res.atoms )

        for aa in residue_atoms:
            coord = aa.coord

        if needs_cterm_capping(res, stru, residue_list):
            residue_atoms.extend(cap_residue(res, stru, cterm_cap, 'cterm'))

        atoms.extend( residue_atoms )


    return StructureRegion( atoms )
