from typing import List

#TODO(CJ): documentation

from copy import deepcopy

import xml.etree.ElementTree as ET

from ..structure import Structure, Solvent, Chain, Residue, Atom

from .api import (
    StructureConstraint,
    CartesianFreeze,
    DistanceConstraint,
    AngleConstraint,
    DihedralConstraint,
    ResiduePairConstraint,
    create_residue_pair_constraint
    )


#TODOC(CJ): generally need more checks and everying in this thing


def structure_constraints_from_xml(topology: Structure, file:str) -> List[StructureConstraint]:
    """TODO(CJ)"""

    tree = ET.parse(file)
    root = tree.getroot()
   

    result = list()
    for rr in root:
        if rr.tag == 'ResiduePairConstraint':
            result.append(residue_pair_constraint_from_xml(topology, rr))
    #TODO(CJ): some kind of check with file system 

    return result


def residue_pair_constraint_from_xml(topology, node):

   
    r1_key = None
    r2_key = None
    r1_atoms = None
    r2_atoms = None
    distanceAB = None
    angle_A = None
    angle_B = None
    torsion_A = None
    torsion_B = None
    torsionAB = None



    for child in node:
        data = child.attrib
        if child.tag == 'Residue1':
            r1_key = (data['chain'], int(data['idx']))
            r1_atoms = tuple(data['atoms'].split(','))

        if child.tag == 'Residue2':
            r2_key = (data['chain'], int(data['idx']))
            r2_atoms = tuple(data['atoms'].split(','))

        if child.tag == 'distanceAB':
            distanceAB = deepcopy(data)

        if child.tag == 'angle_A':
            angle_A = deepcopy(data)

        if child.tag == 'angle_B':
            angle_B = deepcopy(data)

        if child.tag == 'torsion_A':
            torsion_A = deepcopy(data)
        
        if child.tag == 'torsion_B':
            torsion_B = deepcopy(data)

        if child.tag == 'torsionAB':
            torsionAB = deepcopy(data)


    return create_residue_pair_constraint( 
        topology,
        r1_key,
        r2_key,
        r1_atoms,
        r2_atoms,
        distanceAB, 
        angle_A,
        angle_B,
        torsion_A,
        torsion_B,
        torsionAB
    )        
                            

