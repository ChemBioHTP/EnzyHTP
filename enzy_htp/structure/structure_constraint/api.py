"""Defines the StructureConstraint abstract base class as well as the inheriting classes 


This class is one of the data structure of EnzyHTP.
StructureConstraint stands for a coordinate constrain for the structure, including
freeze coordinate and geometry constrain (distance, angle, dihedral).

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-10-28
"""
import copy
from copy import deepcopy
from typing import List, Tuple, Dict, Any

import numpy as np

from enzy_htp.core.logger import _LOGGER
from ..structure import Structure, Solvent, Chain, Residue, Atom

from abc import ABC, abstractmethod


class StructureConstraint(ABC):
    """Abstract base class defining the 

    Attributes:
        atoms_: A List[Atom] containing all the atoms in the constraint.
        target_value_ : The float() target value for the constraint.
        params_: A Dict[str,float] with arbitrary parameters.
    """

    def __init__(self,
                atoms:List[Atom],
                target_value:float,
                params:Dict[str, Any]
                ):
        """Ctor that exits if the supplied number of atoms is NOT correct."""
        
        self.atoms_ = atoms 
        self.target_value_ = target_value 
        self.params_ = deepcopy( params ) 

        self['calc_method'] = params.get('calc_method', 'rosetta')

        if not self.correct_num_atoms():
            _LOGGER.error("Incorrect number of atoms supplied! Exiting...")
            exit( 1 )


    def is_cartesian_freeze(self) -> bool:
        """Is this a cartesian freeze constraint?"""
        return False

    def is_distance_constraint(self) -> bool:
        """Is this a distance constraint?"""
        return False

    def is_angle_constraint(self) -> bool:
        """Is this an angle constraint?"""
        return False

    def is_dihedral_constraint(self) -> bool:
        """Is this a dihedral constraint?"""
        return False

    def is_residue_pair_constraint(self) -> bool:
        """Is this a residue pair constraint?"""
        return False

    @abstractmethod
    def correct_num_atoms(self) -> bool:
        """Does the given constraint have the correct number of atoms?"""
        pass

    def score_energy(self) -> float:
        """
        """
        if self['calc_method'] == 'rosetta':
            penalty:float = self['penalty'] 
            tolerance:float = self['tolerance']
            difference:float = abs(self.current_geometry() - self.target_value)

            if difference <= tolerance:
                return 0.0
            else:
                return penalty * (difference - tolerance)

        else:
            raise TypeError()

        pass

    def change_topology(self, new_topolgy:Structure) -> None:
        assert False
        pass

    @abstractmethod
    def current_geometry(self) -> float:
        """The current geometry of the constrained atoms."""
        pass

    @property
    def atoms(self) -> List[Atom]:
        """Accessor for the atoms involved in the constrained geometry."""
        return self.atoms_

    @atoms.setter
    def atoms(self, val_in : List[Atom] ) -> None:
        """Setter for the atoms in the constrained geometry. Checks for the correct number of atoms."""
        self.atoms_ = val_in
        self.correct_num_atoms()


    def atom_indices(self) -> List[int]: #TODO(CJ): make this an array
        pass

    @property
    def target_value(self) -> float:
        return self.target_value_

    @target_value.setter
    def target_value(self, val_in : float ) -> None:
        self.target_value_ = val_in

    @property
    def params(self) -> Dict[str, Any]:
        return self.params_

    def __getitem__(self, key:str ) -> float:

        if key not in self.params_:
            raise TypeError()
            #TODO(CJ): error here
            pass

        return self.params_[key]        

    def __setitem__(self, key:str, value:Any ) -> None:
        
        self.params_[key] = value 

class CartesianFreeze(StructureConstraint):
    
    def __init__(self, atoms:List[Atom] ):
        super().__init__(self, atoms, 0.0, dict())
     

    def correct_num_atoms(self) -> bool:
        return len(atoms)

    def score_energy(self) -> float:
        return 0.0

    def current_geometry(self) -> float:
        return 0.0


class DistanceConstraint(StructureConstraint):
  
    def correct_num_atoms(self) -> bool:
        return len(self.atoms) == 2

    def current_geometry(self) -> float:
        return self.atoms[0].distance_to(self.atoms[1])


class AngleConstraint(StructureConstraint):


    def correct_num_atoms(self) -> bool:
        return len(self.atoms) == 3


    def current_geometry(self) -> float:
        a = np.array(self.atoms[0].coord)
        b = np.array(self.atoms[1].coord)
        c = np.array(self.atoms[2].coord)
        
        ba = a - b
        bc = c - b
        
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        return np.arccos(cosine_angle)
        

class DihedralConstraint(StructureConstraint):
    

    def correct_num_atoms(self) -> bool:
        return len(self.atoms) == 4

    def current_geometry(self) -> float:
        p0 = np.array(self.atoms[0])
        p1 = np.array(self.atoms[1])
        p2 = np.array(self.atoms[2])
        p3 = np.array(self.atoms[3])

        b0 = -1.0*(p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2

        # normalize b1 so that it does not influence magnitude of vector
        # rejections that come next
        b1 /= np.linalg.norm(b1)

        # vector rejections
        # v = projection of b0 onto plane perpendicular to b1
        #   = b0 minus component that aligns with b1
        # w = projection of b2 onto plane perpendicular to b1
        #   = b2 minus component that aligns with b1
        v = b0 - np.dot(b0, b1)*b1
        w = b2 - np.dot(b2, b1)*b1

        # angle between v and w in a plane is the torsion angle
        # v and w may not be normalized but that's fine since tan is y/x
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        return np.degrees(np.arctan2(y, x))


class ResiduePairConstraint(StructureConstraint):
   

    def __init__(self, 
        residue1:Residue,
        residue2:Residue,
        residue1_atoms:List[Atom],
        residue2_atoms:List[Atom],
        distanceAB:DistanceConstraint=None,
        angle_A:AngleConstraint=None,
        angle_B:AngleConstraint=None,
        torsion_A:DihedralConstraint=None,
        torsion_B:DihedralConstraint=None,
        torsionAB:DihedralConstraint=None
        ):
        
        self.residue1_ = residue1
        self.residue2_ = residue2
        self.residue1_atoms_ = residue1_atoms
        self.residue2_atoms_ = residue2_atoms
        self.distanceAB_ = distanceAB
        self.angle_A_ = angle_A 
        self.angle_B_ = angle_B
        self.torsion_A_ = torsion_A
        self.torsion_B_ = torsion_B
        self.torsionAB_ = torsionAB

        atoms = set()
        for cst in [self.distanceAB_,self.angle_A_, self.angle_B_, self.torsion_A_, self.torsion_B_, self.torsionAB_]:
            #TODO(CJ): should update the 
            if cst is None:
                continue
            #print(cst)
            #print(cst.atoms)
            atoms.update( cst.atoms)
        
        self.atoms_ = atoms
        self.correct_num_atoms()

        #TODO(CJ): should probably assemble atoms and call the correct_num_atoms?
    def is_residue_pair_constraint(self) -> bool:
        """Always True for this class."""
        return True 

    def correct_num_atoms(self) -> bool:
        return len(self.atoms) == 6 

    def current_geometry(self) -> Dict[str, float]:
    
        result:Dict[str,float] = dict()
        
        result['distanceAB'] = self.distanceAB_.current_geometry()
        result['angle_A'] = self.angle_A_.current_geometry()
        result['angle_B'] = self.angle_B_.current_geometry()
        result['torsion_A'] = self.torsion_A_.current_geometry()
        result['torsion_B'] = self.torsion_B_.current_geometry()
        result['torsion_AB'] = self.torsion_AB_.current_geometry()

        return result


    def score_energy(self) -> float:

        energy:float=0.0

        for cst in self.child_constraints:
            energy += cst.score_energy()

        return energy


    @property
    def residue1_atoms(self):
        return self.residue1_atoms_

    @property
    def residue2_atoms(self):
        return self.residue2_atoms_


    @property
    def residue1(self) -> Residue:
        return self.residue1_

    @property
    def residue2(self) -> Residue:
        return self.residue2_

    @property
    def child_constraints(self) -> List[StructureConstraint]:
        return list(filter(
            lambda pr: pr[-1] is not None,
            [('distanceAB', self.distanceAB_),
                ('angle_A', self.angle_A_),
                ('angle_B', self.angle_B_),
                ('torsion_A', self.torsion_A_),
                ('torsion_B', self.torsion_B_),
                ('torsionAB', self.torsionAB_)]
        ))


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
    """
    
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

    residue1:Residue = topology.get_residue(r1_key_str)
    residue2:Residue = topology.get_residue(r2_key_str)

    temp = list()
    for r1 in r1_atoms:
        r1_key = r1_key_str + '.' + r1        
        temp.append(topology.get_atom(r1_key) )

    r1_atoms = temp

    temp = list()
    for r2 in r2_atoms:
        r2_key = r2_key_str + '.' + r2        
        temp.append(topology.get_atom(r2_key) )

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
        
        csts[cidx] = ctor( atom_list, target_value, cst )

    return ResiduePairConstraint(
                                residue1,
                                residue2,
                                r1_atoms,
                                r2_atoms,
                                *csts
                                )





def build_from_preset(topology: Structure,
                      keyword: str,) -> StructureConstraint:
    """constructor that allows building a StructureConstraint from a keyword.
    recommand combining this with functools.partial to make general constrains
    that can apply to different structures
    Args:
        topology: the topology defined by a Structure()
        keyword: the preset constrain keyword. Supported list:
            freeze_backbone: freeze the coordinate of all backbone atoms"""
    supported_keywords = ["freeze_backbone", "freeze_hydrogens"]
    if keyword == "freeze_backbone":
        freeze_atoms = topology.backbone_atoms()
        result = StructureConstraint(topology, freeze_atoms, geom_constrain=[])
        result._is_backbone_freeze = True
    else:
        _LOGGER.error(f"using unsupported keyword {keyword}. (Supported: {supported_keywords})")
        raise ValueError

    return result

