"""Implements a StrutureCluster class which is meant to be essentially a List[Structure] with some additional information about RMSD 
differences between each member in the List. This file contains the StructureCluster class definition as well as the cluster_structures() 
method, a free function which creates StructureCluster()'s from a raw List[Structure].

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2025-07-05
"""
from typing import List

import numpy as np

from .structure import Structure
from enzy_htp import interface


class StructureCluster:
    """Represents a geometric cluster of Structure's that with low RMSDs between constituent members. Functionally, 
    it is a List[Structure] with more information about member-member RMSDs stored in an RMSD matrix.

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

    #TODO(CJ): add structure in method
    def add_structure(self, new_stru:Structure) -> None:
        self.structures_.append( new_stru )

    def structures( self ) -> Structure:
        return self.structures_

    @property
    def rmsd_matrix( self ):
        return self.rmsd_matrix_
    
    def average_rmsd(self, new_stru ):
        values = list()        
        for ss in self.structures():
            values.append( self.rmsd_matrix[new_stru.data['cluster_idx'], ss.data['cluster_idx']]) 
        
        return np.mean(np.array( values ))


def cluster_structures(
    structures:List[Structure],
    align_sele:str,         
    rmsd_sele:str,
    rmsd_dist:float) -> List[StructureCluster]:
    """
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
            if cluster.average_rmsd( ss ) <= rmsd_dist:
                cluster.add_structure( ss )
                break
        else:
            result.append( 
                StructureCluster( rmsd_matrix, ss )
            )
    
    return result
