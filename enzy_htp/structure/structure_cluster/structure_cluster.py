"""Implements a StrutureCluster class which is meant to be essentially a List[Structure] with some additional information about RMSD 
differences between each member in the List. This file contains the StructureCluster class definition as well as the cluster_structures() 
method, a free function which creates StructureCluster()'s from a raw List[Structure].

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2025-07-05
"""

from typing import List

import numpy as np

from enzy_htp.structure import Structure
from enzy_htp import interface


class StructureCluster:
    """Represents a geometric cluster of Structure's that with low RMSDs between constituent members. Functionally, 
    it is a List[Structure] with more information about member-member RMSDs stored in an RMSD matrix. Each Structure is expected to 
    have energeic information for the purpose of ranking StructureCluster's as well as individual member Structure's.

    Attributes:
        structures_: The List[Structure] containing members.
        rmsd_matrix_: The List[List[float]] matrix with RMSD's between all members. Has dimensions of len(structures_)xlen(structures_)
    """
    def __init__(self, rmsd_matrix:List[List[float]], first_stru:Structure=None ):
        """Simple constructor which takes the RMSD matrix and optionally the first Structure to add."""
        self.structures_ = list()
        self.rmsd_matrix_ = rmsd_matrix

        if first_stru:
            self.structures_.append( first_stru )

    def add_structure(self, new_stru:Structure) -> None:
        """Adds a single Structure to the StructureCluster."""
        self.structures_.append( new_stru )

    def structures( self ) -> Structure:    
        """Getter for the member Structure's in the StructureCluster."""
        return self.structures_

    @property
    def rmsd_matrix( self ) -> List[List[float]]:
        """Getter for the RMSD matrix for all Structure's."""
        return self.rmsd_matrix_
    
    def average_rmsd(self, new_stru:Structure ) -> float:
        """Calculates the average RMSD between a candidate Structure and all member Structure's in the 
        StructureCluster object.

        Args:
            new_stru: The candidate Structure.
        
        Returns:
            The averaged RMSD between the candidate Structure and all member Structures.
        """
        values = list()        
        for ss in self.structures():
            values.append( self.rmsd_matrix[new_stru.data['cluster_idx'], ss.data['cluster_idx']]) 
        
        return np.mean(np.array( values ))

   
    def average_score(self) -> float:
        """Gets the average score of the member Structure's"""
        values:List[float] = list()
        for stru in self.structures():
            values.append( stru.data['rosetta_score'] )

        return np.mean(np.array( values ))

    def lowest_energy_structure(self) -> Structure:
        """Which member Structure has the lowest scored energy?"""
        return sorted(
            self.structures(), key=lambda stru: stru.data['rosetta_score']
        )[0]


def cluster_structures(
    structures:List[Structure],
    align_sele:str,         
    rmsd_sele:str,
    rmsd_dist:float) -> List[StructureCluster]:
    """Function that takes a List[Structure] and transforms it into a List[StructureCluster]. Canonical method for 
    creating StructureCluster objects using RMSD as the primary metric.

    Args:
        structures: The raw List[Structure] to be clustered by RMSD.
        align_sele: The pymol-formatted selection to align the constituent Structure's.
        rmsd_sele: The pymol-formatted selection over which Structure-Structure RMSD's are calculated.
        rmsd_dist: The RMSD cutoff at which point a new StructureCluster is created.

    Returns:
        The List[StructureCluster] objects.
    """
    rmsd_matrix = interface.pymol.rmsd_matrix( 
        structures,
        align_sele,
        rmsd_sele
    )
    result:List[StructureCluster] = [ StructureCluster(
        rmsd_matrix, structures.pop() ) ]

    for ss in structures:
        for cluster in result:
            if cluster.average_rmsd(ss) <= rmsd_dist:
                cluster.add_structure( ss )
                break
        else:
            result.append( StructureCluster( rmsd_matrix, ss ) )            
        
               
    return result
