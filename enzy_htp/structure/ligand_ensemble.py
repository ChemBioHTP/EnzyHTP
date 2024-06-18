#TODO(CJ): documentation
from typing import List, Tuple, Dict

from plum import dispatch
from .ligand import (
    Ligand
)

from copy import deepcopy

class LigandEnsemble:


    def __init__(self,
        ligand_0:Ligand
    ):
        self.ligand_0_:Ligand = ligand_0 
        self.ligands_:List[Ligand] = [ligand_0]

    def n_conformers(self) -> int:
        return len(self.ligands())
    
    def ligands(self) -> List[Ligand]:
        return self.ligands_
   
    @property
    def ligand_0(self) -> Ligand:
        return self.ligand_0_

    @dispatch
    def add_conformer(self, new_ligand:Ligand) -> None:

        assert len(new_ligand.atoms) == len(self.ligand_0.atoms)

        new_coords:List[Tuple[float, float, float]] = list()
        for na, oa in zip(new_ligand.atoms, self.ligand_0.atoms):
            assert na.element == oa.element            
            new_coords.append( na.coord )
        
        self.add_conformer( new_coords )

    @dispatch
    def add_conformer(self, coords:List[Tuple[float, float, float]]) -> None:
        
        new_ligand = deepcopy( self.ligand_0 )

        for coord, atom in zip(coords, self.ligand_0.atoms):
            atom.coord = coord

        self.ligands_.append( new_ligand )

    def fix_atom_names(self, name_mapper:Dict[str, str]) -> None:
        for lig in self.ligands_:
            for aa in lig.atoms:
                mapped_name = name_mapper.get( aa.name.strip(), None)
                if mapped_name:
                    aa.name = mapped_name
