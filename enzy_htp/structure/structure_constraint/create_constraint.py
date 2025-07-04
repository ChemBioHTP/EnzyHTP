"""Sub-module containing the constructor functions for various types of StructureConstraint. All of the available
functions are listed below:

    + create_cartesian_freeze()
    + create_backbone_freeze()
    + create_distance_constraint()
    + create_angle_constraint()
    + create_dihedral_constraint()
    + create_hydrogen_bond_freeze()

    + merge_cartesian_freeze()
In practice these functions should be used instead of calling each StructureConstraint()'s constructor
when one wants to make a specific StructureConstraint.

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
2024-01-01
"""

import copy
from copy import deepcopy
from typing import List, Tuple, Dict, Any, Set, Union

import numpy as np

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import WrongTopology
from ..structure import Structure, Solvent, Chain, Residue, Atom
from ..structure_region import StructureRegion, ResidueCap
from ..structure_selection_class import StruSelection
from enzy_htp import config as eh_config

from abc import ABC, abstractmethod


from .api import (
    StructureConstraint,
    CartesianFreeze,
    DistanceConstraint,
    AngleConstraint,
    DihedralConstraint,
    BackBoneFreeze,
    GroupDistanceConstraint,
    )

# region == constructors ==

def create_cartesian_freeze(
        atoms: Union[List[Atom], str],
        topology: Structure= None,) -> CartesianFreeze:
    """Constructor for CartesianFreeze. Takes a List[Atom] or str, as well as optionally a topology.
    Performs basic checks. Note that the str input is not currently supported and will throw.
    
    Args:
        atoms: Either a List[Atom] to be frozen or a str.
        topology: Optionally the Structure to base the constraints on.

    Returns:
        The CartesianFreeze constraint which freezes the supplied atoms.        
    """
    if isinstance(atoms, str):
        if not isinstance(topology, Structure):
            _LOGGER.error("key str types of specification is used but reference "
                          "topology (i.e.: structure) is not supplied. "
                          f"(key: {atoms}, topology: {topology})")
            raise TypeError
        # atoms = stru_sele.select_stru(atoms).atoms
        # TODO figure out how to resolve this loop importing issue
        # structure_constraint -> structure_selection -> interface -> structure_constraint
        # probably make 2 interface component 1 high 1 low?
        raise Exception("TODO")
    result = CartesianFreeze(atoms=atoms)
    result.check_consistent_topology()
    return result

#TODO(CJ): add an overload for the ResidueCap's in StructureRegion

def create_backbone_freeze(stru: Structure, params:Dict=None,) -> BackBoneFreeze:
    """Constructor for BackboneFreeze. Constrains only C,CA,N. Optionally takes a params dict()
    that updates the 
    Args:
        stru: The Structure() whose backbone atoms you want to freeze.
        params: The new params to apply to the BackBoneFreeze. Optional.

    Returns:
        A BackBoneFreeze StructureConstraint for the supplied Structure.
    """
    atoms = stru.backbone_atoms()
    if not atoms:
        return atoms
    result = BackBoneFreeze(atoms=atoms)
    if params is not None:
        result.update_params( params )
    return result

def create_distance_constraint(
        atom_1: Union[Atom, str],
        atom_2: Union[Atom, str],
        target_value: float,
        topology: Structure= None,
        params:Dict=None
        ) -> DistanceConstraint:
    """Constructor for DistanceConstraint between atoms.
    Args:
        atom_1, atom_2:
            the 2 atoms of the distance constraint. they can be specified
            either through Atom() objects or keyword str of Structure().get()
            combining with a Structure() (e.g.: A.100.CA)
        target_value:
            the target distance of the constraint
        topology:
            the reference topology if keyword str is used for atom_1, atom_2.
        params:
            The energy constraint parameters to apply to the DistanceConstraint. Optional.

    Returns:
        The DistanceConstraint object described by the supplied variables.
    """
    atom_1 = _dispatch_get_key(atom_1, topology)
    atom_2 = _dispatch_get_key(atom_2, topology)
    result = DistanceConstraint(
        atoms=[atom_1, atom_2],
        target_value=target_value
    )
    result.check_consistent_topology()

    if params is not None:
        result.update_params(params)
    return result

def create_group_distance_constraint(
        atom_grp_1: str,
        atom_grp_2: str,
        target_value: float,
        topology: Structure= None,
        params:Dict=None
        ) -> GroupDistanceConstraint:
    """Constructor for GroupDistanceConstraint between two groups of atoms.
    Args:
        atom_grp_1, atom_grp_2:
            the 2 atom groups of the distance constraint. they are specified
            through a structure selection pattern in PyMol style.
            (https://pymolwiki.org/index.php/Selection_Algebra)
        target_value:
            the target distance of the constraint
        topology:
            the reference topology if keyword str is used for atom_1, atom_2.
        params:
            The energy constraint parameters to apply to the GroupDistanceConstraint. Optional.

    Returns:
        The GroupDistanceConstraint object described by the supplied variables.
    """
    from ..structure_selection import select_stru # TODO see structure region comment
    atom_grp_1 = select_stru(topology, atom_grp_1).atoms
    atom_grp_2 = select_stru(topology, atom_grp_2).atoms
    result = GroupDistanceConstraint(
        grp_1_atoms=atom_grp_1,
        grp_2_atoms=atom_grp_2,
        target_value=target_value,
    )
    result.check_consistent_topology()

    if params is not None:
        result.update_params(params)
    return result

def create_angle_constraint(
        atom_1: Union[Atom, str],
        atom_2: Union[Atom, str],
        atom_3: Union[Atom, str],
        target_value: float,
        topology: Structure= None,
        params:Dict=None
        ) -> AngleConstraint:
    """Constructor for AngleConstraint between atoms.
    Args:
        atom_1, atom_2, atom_3:
            the 3 atoms of the angle constraint. they can be specified
            either through Atom() objects or keyword str of Structure().get()
            combining with a Structure() (e.g.: A.100.CA)
        target_value:
            the target angle of the constraint
        topology:
            the reference topology if keyword str is used for atom_1, atom_2, atom_3.
        params:
            The energy constraint parameters to apply to the AngleConstraint. Optional.

    Returns:
        The AngleConstraint object described by the supplied variables.
    """
    atom_1 = _dispatch_get_key(atom_1, topology)
    atom_2 = _dispatch_get_key(atom_2, topology)
    atom_3 = _dispatch_get_key(atom_3, topology)
    result = AngleConstraint(
        atoms=[atom_1, atom_2, atom_3],
        target_value=target_value
    )
    result.check_consistent_topology()
    if params is not None:
        result.update_params(params)
    return result

def create_dihedral_constraint(
        atom_1: Union[Atom, str],
        atom_2: Union[Atom, str],
        atom_3: Union[Atom, str],
        atom_4: Union[Atom, str],
        target_value: float,
        topology: Structure= None,
        params:Dict=None
        ) -> DihedralConstraint:
    """Constructor for DihedralConstraint between atoms.
    Args:
        atom_1, atom_2, atom_3, atom_4:
            the 4 atoms of the dihedral constraint. they can be specified
            either through Atom() objects or keyword str of Structure().get()
            combining with a Structure() (e.g.: A.100.CA)
        target_value:
            the target dihedral of the constraint
        topology:
            the reference topology if keyword str is used for
            atom_1, atom_2, atom_3, atom_4.
        params:
            The energy constraint parameters to apply to the DihedralConstraint. Optional.

    Returns:
        The DihedralConstraint object described by the supplied variables.
    """
    atom_1 = _dispatch_get_key(atom_1, topology)
    atom_2 = _dispatch_get_key(atom_2, topology)
    atom_3 = _dispatch_get_key(atom_3, topology)
    atom_4 = _dispatch_get_key(atom_4, topology)
    result = DihedralConstraint(
        atoms=[atom_1, atom_2, atom_3, atom_4],
        target_value=target_value
    )
    result.check_consistent_topology()
    if params is not None:
        result.update_params(params)

    return result

def _dispatch_get_key(target: Union[str, Atom, Residue], stru: Structure) -> Union[Atom, Residue]:
    """apply the key with .get() that support a dispatch of using the
    object (e.g.: Atom()) itself"""
    if isinstance(target, str):
        if not isinstance(stru, Structure):
            _LOGGER.error("key str types of specification is used but reference "
                          "topology (i.e.: structure) is not supplied. "
                          f"(key: {target}, topology: {stru})")
            raise TypeError
        target = stru.get(target)
    return target

# endregion

# region == tools that operations on list/dict of StructureConstraint ==
def merge_cartesian_freeze(cart_freeze_list: List[CartesianFreeze]) -> CartesianFreeze:
    """try to merge a list of cartesian freeze into 1 single
    Cartesian Freeze object. Raise error when they can't merge"""

    # the case only backbone_freeze exists
    temp_iter = iter(cart_freeze_list)
    ref_cons = next(temp_iter)
    try:
        while ref_cons.constraint_type == "backbone_freeze":
            ref_cons = next(temp_iter)
    except StopIteration:
        _LOGGER.debug(f"Only found specialized backbone_freeze. Returning just BackboneFreeze")
        return ref_cons

    # check consistency
    result = CartesianFreeze.clone_from(ref_cons)
    for cons in cart_freeze_list:
        if not cons.is_cartesian_freeze:
            _LOGGER.error(f"this function only take List[CartesianFreeze]. Found: {type(cons)}")
            raise TypeError
        if cons.params != ref_cons.params:
            _LOGGER.error(f"Inconsistent params in the target list. {cons.params} -vs- {ref_cons.params} ")
            raise ValueError
        if cons.topology is not ref_cons.topology:
            _LOGGER.error(f"Inconsistent topology in the target list. {cons.topology} -vs- {ref_cons.topology} ")
            raise ValueError
    # merge atoms
        result.atoms.extend(cons.atoms)
    result.atoms = list(set(result.atoms))
    return result
# endregion


def create_hydrogen_bond_freeze(sr:StructureRegion, params:Dict=None ) -> List[DistanceConstraint]:
    """Creates a List[DistanceConstraint] that constrains all hydrogen bonds in the supplied StructureRegion. 
    Primarily used during QM geometry optimization.

    Args:
        sr: The StructureRegion() where the hydrogen bonds are to be frozen.
        params: The energy parameters for the DistanceConstraint()'s. Optional.

    Returns:
        A List[DistanceConstraint] constraining all hydrogen bonds.
    """
    result:List[DistanceConstraint] = list()
    
    for atom in sr.atoms:
        if atom.element != 'H':
            continue
        
        closest_atom:Atom=sr.closest_heavy_atom(atom)
        target_dist:float=atom.distance_to(closest_atom)
        result.append(create_distance_constraint(atom, closest_atom, target_value=target_dist, params=params))
    
    return result
