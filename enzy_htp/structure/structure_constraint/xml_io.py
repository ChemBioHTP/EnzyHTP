"""Portion of the enzy_htp.structure.structure_constraint submodule that is responsible for creating
StructureConstraint objects from .xml files. This is generally the preferred way to create new StructureConstraint objects.
Functions included are:

    + structure_constraints_from_xml()
    + cartesian_freeze_from_xml()
    + backbone_freeze_from_xml()
    + distance_constraint_from_xml()
    + angle_constraint_from_xml()
    + dihedral_constraint_from_xml()
    + residue_pair_constraint_from_xml()

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-01-01
"""

from typing import List, Dict, Any

from copy import deepcopy

import xml.etree.ElementTree as ET

from ..structure import Structure, Solvent, Chain, Residue, Atom

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs

from .api import (
    StructureConstraint,
    CartesianFreeze,
    BackBoneFreeze,
    DistanceConstraint,
    AngleConstraint, 
    DihedralConstraint,
    ResiduePairConstraint,
    )


from .create_constraint import (
    create_residue_pair_constraint,
    create_cartesian_freeze,
    create_backbone_freeze,
    create_distance_constraint,
    create_angle_constraint,
    create_dihedral_constraint,
    merge_cartesian_freeze,
    )

def structure_constraints_from_xml(topology: Structure, file:str) -> List[StructureConstraint]:
    """Creates a List[StructureConstraint] from a .xml file describing the constrained geometry. The 
    constraints are applied/mapped to a supplied Structure object. Basic checks are performed on the 
    .xml file, the syntax of the supplied constraints, and the compatibility of the Constraints and 
    the Structure object.
    
    Args:
        topology: The Structure object where the StructureConstraint objects will be applied.
        file: The .xml file to get the constraints from.

    Returns:
        A List[StructureConstraint] objects described in the .xml file.
    """

    fs.check_not_empty( file )

    tree = ET.parse(file)
    root = tree.getroot()

    result = list()
    for rr in root:
        if rr.tag == 'ResiduePairConstraint':
            result.append(residue_pair_constraint_from_xml(topology, rr))
        elif rr.tag == 'CartesianFreeze':
            result.append(cartesian_freeze_from_xml(topology, rr))
        elif rr.tag == 'DistanceConstraint':
            result.append(distance_constraint_from_xml(topology, rr))
        elif rr.tag == 'AngleConstraint':
            result.append(angle_constraint_from_xml(topology, rr))
        elif rr.tag == 'DihedralConstraint':
            result.append(dihedral_constraint_from_xml(topology, rr))
        elif rr.tag == 'BackBoneFreeze':
            result.append(backbone_freeze_from_xml(topology,rr))
        else:
            _LOGGER.error(f"The supplied constraint {rr.tag} is not supported!")
            raise TypeError()

    return result

def check_valid_residue_elem(elem:ET.Element) -> None:
    """Helper method that checks if a supplied ET.Element has all the attributes to correctly map
    to a Residue()."""
    
    error:bool = False
    msg:str = 'Residue tag missing attributes: '
    for kw in 'chain idx atoms'.split():
        if kw not in elem:
            msg += kw + ','
            error = True

    msg = msg[:-1]
    if error:
        msg = f"Invalid Residue tag! {msg}"
        _LOGGER.error(msg)
        raise TypeError(msg)


def check_valid_child_constraint(elem:ET.Element) -> None:
    """Helper method that checks if a supplied ET.Element has all the attributes to correctly map
    to a child constraint."""
    
    error = False

    msg = 'Child constraint missing attributes: '
    
    for kw in 'penalty tolerance target_value'.split():
        if kw not in elem:
            mgs += kw + ','
            error = True

    msg = msg[:-1]
    if error:
        _LOGGER.error(f"Invalid Child Constraint tag! {msg}")
        raise TypeError(f"Invalid Child Constraint tag! {msg}")


def residue_pair_constraint_from_xml(topology:Structure, elem: ET.Element) -> ResiduePairConstraint:
    """

    Args:
        topology:
        elem:
    
    Returns:

    """
   
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

    for child in elem:
        data = child.attrib
        if child.tag == 'Residue1':
            check_valid_residue_elem(data)
            r1_key = (data['chain'], int(data['idx']))
            r1_atoms = tuple(data['atoms'].split(','))

        if child.tag == 'Residue2':
            check_valid_residue_elem(data)
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
  

def energy_params_from_element(elem:ET.Element) -> Dict:
    """Parses energy parameters for various packages from the supplied Element. For a given
    package and setting, assumes a format of <package>_<param>=<value>. Returns all results 
    in a nested dict() format. For example, if the element contains the following parameter 
    value pairs: rosetta_tolerance="20.0" amber_rk1="100", the result will be:
    { 'rosetta': {'tolerance': 20.0}, {'amber': {'rk1': 100.0}}

    Note that all values are converted to float()'s. Method deletes the Element attribute-value
    pairs that are used.

    Args:
        elem: The Element to extract energy parameters from

    Returns:
        A dict() with the parsed parameters from the suplied Element.
    """
    to_delete:List[str] = list()
    params:Dict = dict()
    for key_name, key_value in elem.attrib.items():
        key_tks:List[str] = key_name.split('_',1)
        
        if len(key_tks) == 1:
            err_msg:str = f"Invalid format in StructureConstraint parameter {key_tks[0]}. Seperate by underscores!"
            raise TypeError(err_msg)
        
        package_name:str=key_tks[0]
        
        if package_name not in params:
            params[package_name] = dict()
        
        param_name = key_tks[1]
        params[package_name][param_name] = float(key_value)
        to_delete.append( key_name ) 

    for td in to_delete:
        elem.attrib.pop(td)

    return params


def constraint_params_from_element(topology:Structure, elem:ET.Element, n_atoms:int) -> Dict:
    """Helper method that extracts generic values from the supplied Element. They can be used
    to generate all StructureConstraint types except CartesianFreeze and BackBoneFreeze. Deletes
    the Element attribute-value pairs that are used.

    Args:
        topology: The Structure to harvest Atom()'s from. 
        elem: The ElementTree element that values will be extracted from.
        n_atoms: The number of atoms to look for. For N atoms, assumes atom keys atom1, atom2, ... atomN.
        
    Returns:
        A Dict always containing the keys 'atoms', 'target_value', and 'params'.
    """
    to_remove:List[str] = list()
    result:Dict = {'atoms':[], 'target_value':0.}
    for aidx in range(1, 1+n_atoms):
        atom_name:str = f"atom{aidx}"
        atom_key:str = elem.attrib.get(atom_name, None)
        if atom_key is None:
            err_msg:str=f"Expected {n_atoms} atoms in StructureConstraint tag. Could not find attribute '{atom_name}'"
            _LOGGER.error(err_msg)
            raise TypeError(err_msg)
        
        result['atoms'].append(topology.get(atom_key))
        to_remove.append( atom_name )

    target_value_raw:str = elem.attrib.get('target_value', None) 
    to_remove.append('target_value')

    if target_value_raw is None:
        err_msg:str=f"Expected target_value attribute in StructureConstraint tag. Could not find attribute!"
        _LOGGER.error(err_msg)
        raise TypeError(err_msg)

    result['target_value'] = float(target_value_raw)

    for tr in to_remove:
        elem.attrib.pop(tr)

    result['params'] = energy_params_from_element(elem)

    return result

def distance_constraint_from_xml(topology:Structure, elem:ET.Element) -> DistanceConstraint:
    """Creates a DistanceConstraint from an Element found in a .xml file. Performs validation
    checks and removes the attributes that are used from the Element.

    Args:
        topology: The Structure() to bind the DistanceConstraint to.
        elem: The Element object from a .xml file.

    Returns:
        A prepared DistanceConstraint.
    """
    element_params:Dict=constraint_params_from_element(topology, elem, 2)
    return create_distance_constraint( 
        element_params['atoms'][0],
        element_params['atoms'][1],
        element_params['target_value'],
        topology=topology,
        params=element_params['params'])

def angle_constraint_from_xml(topology:Structure, elem:ET.Element) -> AngleConstraint:
    """Creates a AngleConstraint from an Element found in a .xml file. Performs validation
    checks and removes the attributes that are used from the Element.

    Args:
        topology: The Structure() to bind the AngleConstraint to.
        elem: The Element object from a .xml file.

    Returns:
        A prepared AngleConstraint.
    """
    element_params:Dict=constraint_params_from_element(topology, elem, 3)
    return create_angle_constraint( 
        element_params['atoms'][0],
        element_params['atoms'][1],
        element_params['atoms'][2],
        element_params['target_value'],
        topology=topology,
        params=element_params['params'])

def dihedral_constraint_from_xml(topology:Structure, elem:ET.Element) -> DihedralConstraint:
    """Creates a DihedralConstraint from an Element found in a .xml file. Performs validation
    checks and removes the attributes that are used from the Element.

    Args:
        topology: The Structure() to bind the DihedralConstraint to.
        elem: The Element object from a .xml file.

    Returns:
        A prepared DihedralConstraint.
    """
    element_params:Dict=constraint_params_from_element(topology, elem, 4)
    return create_dihedral_constraint( 
        element_params['atoms'][0],
        element_params['atoms'][1],
        element_params['atoms'][2],
        element_params['atoms'][3],
        element_params['target_value'],
        topology=topology,
        params=element_params['params'])


def cartesian_freeze_from_xml(topology:Structure, elem:ET.Element) -> CartesianFreeze:
    """Creates a CartesianFreeze from an Element found in a .xml file. Performs validation
    checks and removes the attributes that are used from the Element.

    Args:
        topology: The Structure() to bind the CartesianFreeze to.
        elem: The Element object from a .xml file.

    Returns:
        A prepared CartesianFreeze.
    """
    raw_atom_str:str=elem.attrib.get('atoms', None)
    if raw_atom_str is None:
        err_msg:str="Did not find attribute atoms! Expected in CartesianFreeze Element"
        _LOGGER.error(err_msg)
        raise TypeError(err_msg)
   
    atoms:List[Atom] = list()
    for atom_tk in raw_atom_str.split(','):
        atoms.append(topology.get(atom_tk))

    elem.attrib.pop('atoms')

    cf = CartesianFreeze(atoms)
    cf.update_params( energy_params_from_element(elem) )
    return cf 


def backbone_freeze_from_xml(topology:Structure, elem:ET.Element) -> BackBoneFreeze:
    """Creates a BackBoneFreeze from an Element found in a .xml file. Performs validation
    checks and removes the attributes that are used from the Element.

    Args:
        topology: The Structure() to bind the BackBoneFreeze to.
        elem: The Element object from a .xml file.

    Returns:
        A prepared BackBoneFreeze.
    """
    return create_backbone_freeze(topology, params=energy_params_from_element(elem ))

