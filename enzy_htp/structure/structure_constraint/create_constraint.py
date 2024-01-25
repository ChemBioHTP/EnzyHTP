"""

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
from enzy_htp import config as eh_config

from abc import ABC, abstractmethod


from .api import (
    StructureConstraint,
    CartesianFreeze,
    DistanceConstraint,
    AngleConstraint,
    DihedralConstraint,
    ResiduePairConstraint,
    BackBoneFreeze
    )



# region == constructors ==
def create_residue_pair_constraint(
                topology:Structure,
                r1_key:Tuple[str,int],
                r2_key:Tuple[str,int],
                r1_atoms:Tuple[str,str,str],
                r2_atoms:Tuple[str,str,str],
                distanceAB:Dict[str,float]=None,
                angle_A:Dict[str,float]=None,
                angle_B:Dict[str,float]=None,
                torsion_A:Dict[str,float]=None,
                torsion_B:Dict[str,float]=None,
                torsionAB:Dict[str,float]=None) -> ResiduePairConstraint:
    """Converts specified 
    
    Args:
        topology:
        r1_key:
        r2_key:
        r1_atoms:
        r2_atoms:
        distanceAB:
        angle_A:
        angle_B:
        torsion_A:
        torsion_B:
        torsionAB:

    Returns:
        A ResiduePairConstraint object.
    """
    #TODO(CJ): a bunch of checks
    r1_key_str = ".".join(map(str,r1_key))
    r2_key_str = ".".join(map(str,r2_key))

    residue1:Residue = topology.get(r1_key_str)
    residue2:Residue = topology.get(r2_key_str)

    temp = list()
    for r1 in r1_atoms:
        r1_key = r1_key_str + '.' + r1        
        temp.append(topology.get(r1_key) )

    r1_atoms = temp

    temp = list()
    for r2 in r2_atoms:
        r2_key = r2_key_str + '.' + r2        
        temp.append(topology.get(r2_key) )

    r2_atoms = temp

    csts = [distanceAB, angle_A, angle_B, torsion_A, torsion_B, torsionAB]
    for cidx, cst in enumerate(csts):
        if cst is None:
            continue

        for kw in 'target_value penalty tolerance'.split():
            assert kw in cst #TODO(CJ): put some kind of error code here

        for cst_kk in cst.keys():
            cst[cst_kk] = float(cst[cst_kk])

        ctor:StructureConstraint = None
        atom_list:List[Atom] = None
        target_value:float = float(cst.pop('target_value'))
        if cidx == 0:
            atom_list = [r1_atoms[0], r2_atoms[0]]
            ctor = DistanceConstraint
        elif cidx <= 2:
            if cidx == 1:
                atom_list = [r1_atoms[1], r1_atoms[0], r2_atoms[0]]
            else:
                atom_list = [r1_atoms[0], r2_atoms[0], r2_atoms[1]]
            ctor = AngleConstraint
        elif cidx <= 5:
            if cidx == 3:
                atom_list = [r1_atoms[2], r1_atoms[1], r1_atoms[0], r2_atoms[0]]
            elif cidx == 4:
                atom_list = [r1_atoms[0], r2_atoms[0], r2_atoms[1], r2_atoms[2]]
            else:
                atom_list = [r1_atoms[1], r1_atoms[0], r2_atoms[0], r2_atoms[1]]
            ctor = DihedralConstraint
        
        csts[cidx] = ctor( atom_list, target_value )

    return ResiduePairConstraint(
                                residue1,
                                residue2,
                                r1_atoms,
                                r2_atoms,
                                *csts
                                )

def create_cartesian_freeze(
        atoms: Union[Atom, str],
        topology: Structure= None,) -> CartesianFreeze:
    """constructor for CartesianFreeze."""
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

def create_backbone_freeze(stru: Structure, params:Dict) -> BackBoneFreeze:
    """constructor for BackboneFreeze. Constrain only C,CA,N.
    Args:
        stru: the target structure"""
    atoms = stru.backbone_atoms()
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
    """constructor for DistanceConstraint between atoms.
    Args:
        atom_1, atom_2:
            the 2 atoms of the distance constraint. they can be specified
            either through Atom() objects or keyword str of Structure().get()
            combining with a Structure() (e.g.: A.100.CA)
        target_value:
            the target distance of the constraint
        topology:
            the reference topology if keyword str is used for atom_1, atom_2."""
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

def create_angle_constraint(
        atom_1: Union[Atom, str],
        atom_2: Union[Atom, str],
        atom_3: Union[Atom, str],
        target_value: float,
        topology: Structure= None,
        params:Dict=None
        ) -> AngleConstraint:
    """constructor for AngleConstraint between atoms.
    Args:
        atom_1, atom_2, atom_3:
            the 3 atoms of the angle constraint. they can be specified
            either through Atom() objects or keyword str of Structure().get()
            combining with a Structure() (e.g.: A.100.CA)
        target_value:
            the target angle of the constraint
        topology:
            the reference topology if keyword str is used for atom_1, atom_2, atom_3."""
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
    """constructor for DihedralConstraint between atoms.
    Args:
        atom_1, atom_2, atom_3, atom_4:
            the 4 atoms of the dihedral constraint. they can be specified
            either through Atom() objects or keyword str of Structure().get()
            combining with a Structure() (e.g.: A.100.CA)
        target_value:
            the target dihedral of the constraint
        topology:
            the reference topology if keyword str is used for
            atom_1, atom_2, atom_3, atom_4."""
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
