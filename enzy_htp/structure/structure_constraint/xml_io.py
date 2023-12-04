from typing import List

#TODO(CJ): documentation

from copy import deepcopy

import xml.etree.ElementTree as ET

from ..structure import Structure, Solvent, Chain, Residue, Atom

from enzy_htp.core import _LOGGER

from .api import (
    StructureConstraint,
    CartesianFreeze,
    DistanceConstraint,
    AngleConstraint, DihedralConstraint,
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

def check_valid_residue_node(data) -> None:
    
    error = False
    msg = 'Residue tag missing attributes: '
    for kw in 'chain idx atoms'.split():
        if kw not in data:
            msg += kw + ','
            error = True

    msg = msg[:-1]
    if error:
        _LOGGER.error(f"Invalid Residue tag! {msg}")
        _LOGGER.error("Exiting...")
        exit( 1 )

def check_valid_child_constraint(data) -> None:
    
    error = False

    msg = 'Child constraint missing attributes: '
    
    for kw in 'penalty tolerance target_value'.split():
        if kw not in data:
            mgs += kw + ','
            error = True

    msg = msg[:-1]
    if error:
        _LOGGER.error(f"Invalid Child Constraint tag! {msg}")
        _LOGGER.error("Exiting...")
        exit( 1 )


def residue_pair_constraint_from_xml(topology, node) -> ResiduePairConstraint:

        
   
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
            check_valid_residue_node(data)
            r1_key = (data['chain'], int(data['idx']))
            r1_atoms = tuple(data['atoms'].split(','))

        if child.tag == 'Residue2':
            check_valid_residue_node(data)
            r2_key = (data['chain'], int(data['idx']))
            r2_atoms = tuple(data['atoms'].split(','))

        if child.tag == 'distanceAB':
            check_valid_child_constraint(data)
            distanceAB = deepcopy(data)

        if child.tag == 'angle_A':
            check_valid_child_constraint(data)
            angle_A = deepcopy(data)

        if child.tag == 'angle_B':
            check_valid_child_constraint(data)
            angle_B = deepcopy(data)

        if child.tag == 'torsion_A':
            check_valid_child_constraint(data)
            torsion_A = deepcopy(data)
        
        if child.tag == 'torsion_B':
            check_valid_child_constraint(data)
            torsion_B = deepcopy(data)

        if child.tag == 'torsionAB':
            check_valid_child_constraint(data)
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
                            

