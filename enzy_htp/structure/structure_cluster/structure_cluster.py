
from typing import List

import numpy as np

from enzy_htp.structure import Structure
from enzy_htp import interface


class StructureCluster:


    def __init__(self, rmsd_matrix, first_stru:Structure=None ):
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

   
    #TODO(CJ): add ability to specify scoring
    def average_score(self) -> float:
        
        values:List[float] = list()
        for stru in self.structures():
            values.append( stru.data['rosetta_score'] )

        return np.mean(np.array( values ))

    def lowest_energy_structure(self) -> Structure:

        return sorted(
            self.structures(), key=lambda stru: stru.data['rosetta_score']
        )[0]


def cluster_structures(
    structures:List[Structure],
    align_sele:str,         
    rmsd_sele:str,
    rmsd_dist:float) -> List[StructureCluster]:
    
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
