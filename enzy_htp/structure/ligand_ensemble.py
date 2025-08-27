"""
Defines a LigandEnsemble class that represents a collection of different geometries for the same Ligand. This class manages 
multiple conformers of that single type of ligand, primarily in a manner which is compatible with RosettaLigand and Rosetta 
rotamer generation logic more broadly.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2025-07-05
"""
from typing import List, Tuple, Dict

from plum import dispatch
from .ligand import (
    Ligand
)

from copy import deepcopy

class LigandEnsemble:
    """Represents a collection of Ligand objects, specifically designed to manage conformers of 
    a single ligand type. It allows for the addition of new conformers and provides functionality 
    for manipulating atom names across the ensemble.

    Attributes:
        ligand_0_: The reference Ligand object from which conformers will be derived.
        ligands_: A list containing all geometries in the ensemble.
    """
    def __init__(self,
        ligand_0:Ligand
    ):
        """Simple constructor for LigandEnsemble, takes only an argument for ligand_0."""
        self.ligand_0_:Ligand = ligand_0 
        self.ligands_:List[Ligand] = [ligand_0]

    def n_conformers(self) -> int:
        """How many conformers are present in the ensemble?"""
        return len(self.ligands())
    
    def ligands(self) -> List[Ligand]:
        """Getter for all Ligand conformers."""
        return self.ligands_
   
    @property
    def ligand_0(self) -> Ligand:
        """Getter for the reference Ligand."""
        return self.ligand_0_

    @dispatch
    def add_conformer(self, new_ligand:Ligand) -> None:
        """Adds a conformer to the ensemble using a new Ligand object.

        Args:
            new_ligand: The Ligand object representing the new conformer added to the ensemble.

        Returns:
            Nothing.

        Raises:
            AssertionError: If the new ligand does not have the same number of atoms as the reference ligand.
        """
        assert len(new_ligand.atoms) == len(self.ligand_0.atoms)

        new_coords:List[Tuple[float, float, float]] = list()
        for na, oa in zip(new_ligand.atoms, self.ligand_0.atoms):
            assert na.element == oa.element            
            new_coords.append( na.coord )
        
        self.add_conformer( new_coords )

    @dispatch
    def add_conformer(self, coords:List[Tuple[float, float, float]]) -> None:
        """Adds a new conformer to the ensemble based on coordinates.

        Args:
            coords: A list of tuples representing the new atom coordinates.

        Returns:
            Nothing.

        Generates a deepcopy of the reference ligand and updates its atom coordinates.
        """        
        new_ligand = deepcopy( self.ligand_0 )

        for coord, atom in zip(coords, self.ligand_0.atoms):
            atom.coord = coord

        self.ligands_.append( new_ligand )

    def fix_atom_names(self, name_mapper:Dict[str, str]) -> None:
        """Updates atom names in all ligands based on a provided mapping. This method iterates through each ligand in 
        the ensemble and updates atom names if they match keys in the name_mapper.
        
        Args:
            name_mapper: A dictionary mapping old atom names to new atom names.

        Returns:
            Nothing.
        """
        for lig in self.ligands_:
            for aa in lig.atoms:
                mapped_name = name_mapper.get( aa.name.strip(), None)
                if mapped_name:
                    aa.name = mapped_name
