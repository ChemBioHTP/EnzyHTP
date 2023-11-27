"""Defines the StructureConstraint class.
This class is one of the data structure of EnzyHTP.
StructureConstraint stands for a coordinate constrain for the structure, including
freeze coordinate and geometry constrain (distance, angle, dihedral).

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2022-10-28
"""
import copy
from typing import List, Tuple

from enzy_htp.core.logger import _LOGGER
from ..structure import Structure, Solvent, Chain, Residue, Atom

from abc import ABC


class GeometricConstraint(ABC):

    def __init__(self):
        pass
                                                        

    def score_energy(self) -> float:
        pass


class CartesianFreeze(GeometricConstraint):
    
    def __init__(self, 
     


class AngleConstraint:
    
    def __init__(self, 



class StructureConstraint:
    """This class describe a constrain of a structure (defined by the topology)
    Constructors:
        from_keyword()
    Attributes:
        topology: Structure
        freeze_atoms: List[Atom]
        geom_constrain: List[Tuple[List[Atom], float]]

    Derived properties:
        structures: List[Structure]"""

    def __init__(self, topology: Structure,
                 freeze_atoms: List[Atom], geom_constrain: List[Tuple[List[Atom], float]]):
        self._topology = topology
        self._freeze_atoms = freeze_atoms
        self._geom_constrain = geom_constrain
        # init for check
        self._is_backbone_freeze = None

        #self.distance


    def add_angle_constraint(self, 


    def set_topology(self, topology:Structure) -> None:
        pass

    def is_backbone_freeze(self) -> bool:
        """this method saves effort when parsing to common Amber constrain"""
        if self._is_backbone_freeze is None:
            return set(self._topology.backbone_atoms()) == set(self._freeze_atoms)
        else:
            return self._is_backbone_freeze

    def is_compatible(self, stru: Structure) -> bool:
        """TODO(CJ)"""
        pass


    def energy_penalty(self) -> float:
        pass



def build_from_xml(file:str) -> StructureConstraint:
    pass


def build_from_preset(topology: Structure,
                      keyword: str,) -> StructureConstraint:
    """constructor that allows building a StructureConstraint from a keyword.
    recommand combining this with functools.partial to make general constrains
    that can apply to different structures
    Args:
        topology: the topology defined by a Structure()
        keyword: the preset constrain keyword. Supported list:
            freeze_backbone: freeze the coordinate of all backbone atoms"""
    supported_keywords = ["freeze_backbone"]
    if keyword == "freeze_backbone":
        freeze_atoms = topology.backbone_atoms()
        result = StructureConstraint(topology, freeze_atoms, geom_constrain=[])
        result._is_backbone_freeze = True
    else:
        _LOGGER.error(f"using unsupported keyword {keyword}. (Supported: {supported_keywords})")
        raise ValueError

    return result

