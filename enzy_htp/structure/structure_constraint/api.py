"""Defines the StructureConstraint class.
This class is one of the data structure of EnzyHTP.
StructureConstraint stands for a coordinate constrain for the structure, including
freeze coordinate and geometry constrain (distance, angle, dihedral).

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
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

    def __init__(self,
                atoms:List[Atom],
                target_value:float,
                params:Dict[str, Any]
                ):
        
        self.atoms_ = atoms 
        self.target_value_ = target_value 
        self.params_ = deepcopy( params ) 

        self['calc_method'] = params.get('calc_method', 'rosetta')

        assert self.correct_num_atoms()

    #TODO(CJ): make the is_<type> functions

    @abstractmethod
    def correct_num_atoms(self) -> bool:
        pass

    def score_energy(self) -> float:
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

    def change_topology(self, new_topolgy:Structure)-> None:
        assert False
        pass

    @abstractmethod
    def current_geometry(self) -> float:
        pass

    @property
    def atoms(self) -> List[Atom]:
        return self.atoms_

    @atoms.setter
    def atoms(self, val_in : List[Atom] ) -> None:
        self.atoms_ = val_in


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

    def score_energy(self) -> float:
        if self['calc_method'] == 'rosetta':
            tolerance:float = self['tolerance']
            difference:float = self.deviation()

            if difference <= tolerance:
                return 0.0
            else:
                return penalty * (difference - tolerance)

        else:
            raise TypeError()

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
        


#TODO(CJ): this
class DihedralConstraint(StructureConstraint):
    

    def correct_num_atoms(self) -> bool: return len(self.atoms) == 4

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
        distanceAB:DistanceConstraint=None,
        angle_A:AngleConstraint=None,
        angle_B:AngleConstraint=None,
        torsion_A:DihedralConstraint=None,
        torsion_B:DihedralConstraint=None,
        torsionAB:DihedralConstraint=None
        ):
        
        self.residue1_ = residue1
        self.residue2_ = residue2
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
       
        #print(len(atoms))

        #TODO(CJ): should probably assemble atoms and call the correct_num_atoms?


    def correct_num_atoms(self) -> bool:
        return len(self.atoms) == 6 

    def current_geometry(self) -> Dict[str, float]:
    
        result = dict()
        self.distanceAB_ = distanceAB
        self.angle_A_ = angle_A 
        self.angle_B_ = angle_B
        self.torsion_A_ = torsion_A
        self.torsion_B_ = torsion_B
        self.torsion_AB_ = torsion_AB


    def score_energy(self) -> float:
        assert False


    @property
    def residue1(self) -> Residue:
        return self.residue1_

    @property
    def residue2(self) -> Residue:
        return self.residue2_


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

        ctor:StructureConstraint = None
        atom_list:List[Atom] = None
        target_value:float = cst.pop('target_value')
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

