#TODO(CJ):





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

    @property
    def atoms(self) -> List[Atom]:
        return self.atoms_

    def get_atom(self, key:str ) -> Atom:
        #NOTE(CJ): need this part for the 
        pass


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
