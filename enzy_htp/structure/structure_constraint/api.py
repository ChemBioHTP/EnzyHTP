"""Defines the StructureConstraint abstract base class as well as the inheriting classes that specialize it. 
StructureConstraints are meant to define flexible, non-package specific relationships that can be translated 
in between different software packages. Factory functions for constraints are also included. The contained 
classes are:

    + StructureConstraint: The abstract base class.
    + CartesianFreeze: Specific atoms will be frozen.
    + DistanceConstraint: Constraint that specifies a distance.
    + AngleConstraint: Constraint that specifies an Angle.
    + DihedralConstraint: Constraint that specifies a dihedral.
    + ResiduePairConstraint: Specifies a constraint set between two residues in a style similar to Rosetta EnzDes.

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-10-28
"""
import copy
from copy import deepcopy
from typing import List, Tuple, Dict, Any, Set

import numpy as np

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import WrongTopology
from ..structure import Structure, Solvent, Chain, Residue, Atom
from enzy_htp import config as eh_config

from abc import ABC, abstractmethod


class StructureConstraint(ABC):
    """Abstract base class defining the API for a constraint. Each primitive StructureConstraint defines exactly one 
    type of interaction. StructureConstraints are meant to define flexible, non-package specific relationships that
    can be translated in between different software packages.

    Attributes:
        atoms_:
            A List[Atom] containing all the atoms in the constraint.
        target_value_ :
            The float() target value for the constraint.
        params_:
            A Dict[str,float] with a that contain parameters for EVERY software
            that this constraint type support. user can just edit this mapper to
            change the parameter.
            example:
            {
            "rosetta" : {"XXX" : ###},
            "amber" : {"XXX" : ###},
            "xtb" : {"XXX" : ###},
            }
    """

    def __init__(self,
                 atoms: List[Atom],
                 target_value: float,
                 params: Dict[str, Dict[str, Any]],
                 ):
        """Ctor that exits if the supplied number of atoms is NOT correct."""
        
        self.atoms_ = atoms 
        self.target_value_ = target_value 
        self.params_ = deepcopy( params ) 

        if not self.correct_num_atoms():
            _LOGGER.error("Incorrect number of atoms supplied! Exiting...")
            exit( 1 )

    # region == property getter ==
    @property
    def atoms(self) -> List[Atom]:
        """Accessor for the atoms involved in the constrained geometry."""
        return self.atoms_

    @atoms.setter
    def atoms(self, val_in : List[Atom] ) -> None:
        """Setter for the atoms in the constrained geometry. Checks for the correct number of atoms, errors if it is the incorrect number.."""
        self.atoms_ = val_in
        if not self.correct_num_atoms():
            _LOGGER.error(f"An incorrect number of atoms ({len(self.atoms)})was supplied! Exiting...")
            exit( 1 )

    @property
    def target_value(self) -> float:
        """Getter for the target value of the constrained geometry."""
        return self.target_value_

    @target_value.setter
    def target_value(self, val_in : float ) -> None:
        """Change the stored target value of the constrained geometry.."""
        self.target_value_ = val_in

    @property
    def params(self) -> Dict[str, Any]:
        """the mapper that contain parameters for all support software
        of this type of constraint"""
        return self.params_

    @property
    def topology(self) -> Structure:
        """get the current topology-context of the constraint"""
        self.check_consistent_topology()
        return self.atoms[0].root

    @property
    def constraint_type(self) -> str:
        """hard coded constraint type"""
        return "general"

    def clone(self) -> "StructureConstraint":
        """Clones the StructureConstraint with the original target_value."""
        return self.__init__(self.atoms, self.target_value, deepcopy(self.params_))

    def clone_current(self) -> "StructureConstraint":
        """Get a version of the StructureConstraint with the current value as the 
        target value."""
        return self.__init__(self.atoms, self.current_value, deepcopy(self.params_))

    @abstractmethod
    def current_geometry(self) -> float:
        """The current geometry of the constrained atoms."""
        pass

    def atom_indices(self) -> np.ndarray: 
        """Get the indices of the atoms in the constraint.""" 
        return np.array(
            [aa.idx for aa in self.atoms]
        )
    
    @property
    def atom_names(self) -> Set[str]:
        """get all unique atom names of self.atoms"""
        result = set()
        for atom in self.atoms:
            result.add(atom.name)
        return result
    # endregion

    # region == checker ==
    def check_consistent_topology(self) -> None:
        """check whether Atom()s are from the same Structure()
        raise an error if not."""
        top = self.atoms[0].root
        for atom in self.atoms:
            current_top = atom.root
            if not isinstance(current_top, Structure):
                _LOGGER.error(
                    f"Topology should be Structure(). ({atom} has {current_top})")
                raise WrongTopology

            if current_top is not top:
                _LOGGER.error(
                    "Atom()s in StructureConstraint() "
                    "have inconsistent topology! "
                    f"({self.atoms[0]} has {top} --vs-- {atom} has {current_top})")
                raise WrongTopology

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
    # endregion

    # region == editor ==
    def change_topology(self, new_topology:Structure) -> None:
        """Switch the associated topology and atoms for a given Constraint. Maps the atom keys
        as described in Atom.key() and Structure.get().

        Args:
            new_topology: new structure to translate atoms to.

        Returns:
            Nothing.
        """
        new_atoms:List[Atom] = list()
        for atom in self.atoms:
            target_key:str=atom.key
            new_atoms.append( new_topology.get(target_key))
    
        self.atoms = new_atoms

    # endregion

    # region == special ==
    def __getitem__(self, key:str ) -> float:
        """Bracket getter for the .params dict() in the structure constraint. Exits if the 
        supplied key is not found in the .params dict().
        
        Args:
            key: The str() key for the attribute we are trying to access.

        Returns:
            Value corresponding to the specified key.
        """

        if key not in self.params_:
            _LOGGER.error(f"The attribute '{key}' is not in the constrained geometry params. Exiting...")
            exit( 1 )
        
        return self.params_[key]        

    def __setitem__(self, key:str, value:Any ) -> None:
        """Bracket setter for the .params dict() in the structure constraint. Does not check if the supplied 
        key exists and will overwrite existing values. 

        Args:
            key: The str() key to add to the .params dict().
            value: Value to set in the .params dict().
        
        Returns:
            Nothing.
        """
        self.params_[key] = value 
    # endregion

    #TODO(CJ): add function that checks if topology and constraints are compatible
    #TODO(CJ): will need to make a version of this that actually works for the ResiduePairConstraint
    
    # TODO move this to interface e.g.: RosettaInterface for current
    def score_energy(self) -> float: 
        """Scores the energy penalty of the constraint. If the calc method is not specified, the method exits."""
        if self['calc_method'] == 'rosetta':
            penalty:float = self['penalty'] 
            tolerance:float = self['tolerance']
            difference:float = abs(self.current_geometry() - self.target_value)

            if difference <= tolerance:
                return 0.0
            else:
                return penalty * (difference - tolerance)

        else:
            _LOGGER.error(f"Unspecified calculation method {self['calc_method']}. Supported include: rosetta. Exiting...")
            exit( 1 )


class CartesianFreeze(StructureConstraint):
    """Specialization of StructureConstraint() for Atoms() that are frozen in Cartesian space. Many
    of the methods yield junk values since no measurement is needed for this Class.
    """

    DEFAULT_PARAMS = {
        "amber": {
            "restraint_wt": eh_config["amber.DEFAULT_CARTESIAN_FREEZE_WEIGHT"],
        },
    }
    
    def __init__(self, atoms:List[Atom]):
        super().__init__(atoms, 0.0, self.DEFAULT_PARAMS)

    @property
    def constraint_type(self) -> str:
        """hard coded constraint type"""
        return "cartesian_freeze"

    def is_cartesian_freeze(self) -> bool:
        """Is this a cartesian freeze constraint? Always True for this class."""
        return True

    def correct_num_atoms(self) -> bool:
        """True as long as there is at least one atom in the constraint."""
        return len(self.atoms)

    def current_geometry(self) -> float:
        """Always returns 0.0. Not needed for this type of constraint.
        TODO see if there will be a need"""
        return 0.0


class DistanceConstraint(StructureConstraint):
    """Specialization of StructureConstraint() for the distance between two Atom()'s."""

    @property
    def constraint_type(self) -> str:
        """hard coded constraint type"""
        return "distance_constraint"

    def is_distance_constraint(self) -> bool:
        """Is this a distance constraint? Always True for this class."""
        return True

    def correct_num_atoms(self) -> bool:
        """A distance constraint should have two atoms."""
        return len(self.atoms) == 2

    def current_geometry(self) -> float:
        """The cartesian distance between two 3D points."""
        return self.atoms[0].distance_to(self.atoms[1])


class AngleConstraint(StructureConstraint):
    """Specialization of StructureConstraint() for the angle between three Atom()'s."""

    @property
    def constraint_type(self) -> str:
        """hard coded constraint type"""
        return "angle_constraint"

    def is_angle_constraint(self) -> bool:
        """Is this an angle constraint? Always True for this class."""
        return True 

    def correct_num_atoms(self) -> bool:
        """An angle constraint should have three atoms."""
        return len(self.atoms) == 3

    def current_geometry(self) -> float:
        """The angle between three 3D points using dot-product."""
        return self.atoms[0].angle_with(self.atoms[1], self.atoms[2])
        

class DihedralConstraint(StructureConstraint):
    """Specialization of StructureConstraint() for the dihedral angle between four Atom()'s"""

    @property
    def constraint_type(self) -> str:
        """hard coded constraint type"""
        return "dihedral_constraint"

    def is_dihedral_constraint(self) -> bool:
        """Is this a dihedral constraint?"""
        return True 

    def correct_num_atoms(self) -> bool:
        """A dihedral constraint should have 4 atoms."""
        return len(self.atoms) == 4

    def current_geometry(self) -> float:
        """Measurement for the current dihedral between the 4 atoms in the constraint."""
        return self.atoms[0].dihedral_with(self.atoms[1], self.atoms[2], self.atoms[3])


class ResiduePairConstraint(StructureConstraint):
    """Specialization of StructureConstraint() representing 

    Attributes:
        
    """
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
        for (cst_name,cst) in self.child_constraints:
            if cst is None:
                continue
            atoms.update( cst.atoms)
        
        self.atoms_ = list(atoms)
        self.correct_num_atoms()

    @property
    def constraint_type(self) -> str:
        """hard coded constraint type"""
        return "residue_pair_constraint"


    def clone(self) -> "ResiduePairConstraint":
        """TODO(CJ)"""
        return ResiduePairConstraint.__init__(
                self.residue1_, 
                self.residue2_,
                self.residue1_atoms_,
                self.residue2_atoms_,
                self.distanceAB_,
                self.angle_A_,
                self.angle_B_,
                self.torsion_A_,
                self.torsion_B_,
                self.torsionAB_
                )
        

    def change_topology(self, new_topology: Structure) -> None:
        """Same as StructureConstraint.change_topology() but with extra steps because of the composite 
        nature of the class. First, the residue and residue atoms are mapped over. Next, the child constraints
        are mapped over. Last, the atoms are mapped over."""
        
        self.residue1_ = new_topology.get(self.residue1_.key_str)
        self.residue2_ = new_topology.get(self.residue2_.key_str)
        self.residue1_atoms_ = [new_topology.get(aa.key) for aa in self.residue1_atoms_]
        self.residue2_atoms_ = [new_topology.get(aa.key) for aa in self.residue2_atoms_]

        for (cst_name, cst) in self.child_constraints:
            if cst is not None:
                cst.change_topology(new_topology)

        atoms = set()
        for (cst_name,cst) in self.child_constraints:
            if cst is None:
                continue
            atoms.update( cst.atoms)
        
        self.atoms_ = list(atoms)
        self.correct_num_atoms()

        #TODO(CJ): should probably assemble atoms and call the correct_num_atoms?
    def is_residue_pair_constraint(self) -> bool:
        """Always True for this class."""
        return True 

    def correct_num_atoms(self) -> bool:
        """This type of composite residue needs 6 total atoms."""
        return len(self.atoms) == 6 

    def current_geometry(self) -> Dict[str, float]:
        """Current geometry for all of the child constraints. Packages the results in a dict() of 
        (key, value) format (cst_name, value), where cst_name is a str with Rosetta EnzDes format.
        These include: distanceAB, angle_A, angle_B, torsion_A, torsion_B, torsion_AB.
        """
    
        result:Dict[str,float] = dict()
        
        result['distanceAB'] = self.distanceAB_.current_geometry()
        result['angle_A'] = self.angle_A_.current_geometry()
        result['angle_B'] = self.angle_B_.current_geometry()
        result['torsion_A'] = self.torsion_A_.current_geometry()
        result['torsion_B'] = self.torsion_B_.current_geometry()
        result['torsion_AB'] = self.torsion_AB_.current_geometry()

        return result


    def score_energy(self) -> float:
        """Composite score energy for the """

        energy:float=0.0

        for (cst_name,cst) in self.child_constraints:
            energy += cst.score_energy()

        return energy

    def clone_current(self) -> "ResiduePairConstraint":
        """TODO(CJ)"""
        assert False

    @property
    def residue1_atoms(self):
        """Getter for the 3 Atom() objects in residue 1."""
        return self.residue1_atoms_

    @property
    def residue2_atoms(self):
        """Getter for the 3 Atom() objects in residue 2."""
        return self.residue2_atoms_

    @property
    def residue1(self) -> Residue:
        """Getter for residue 1."""
        return self.residue1_

    @property
    def residue2(self) -> Residue:
        """Getter for residue 2."""
        return self.residue2_

    @property
    def child_constraints(self) -> List[StructureConstraint]:
        """Gets the child contraints that are not None.TODO(CJ) 

        """
        return list(filter(
            lambda pr: pr[-1] is not None,
            [('distanceAB', self.distanceAB_),
                ('angle_A', self.angle_A_),
                ('angle_B', self.angle_B_),
                ('torsion_A', self.torsion_A_),
                ('torsion_B', self.torsion_B_),
                ('torsionAB', self.torsionAB_)]
        ))


class BackBoneFreeze(CartesianFreeze):
    """Specialization of CartesianFreeze() for constraining the movement of
    backbone Atom()s"""

    DEFAULT_PARAMS = {
        "amber": {
            "restraint_wt": eh_config["amber.DEFAULT_CARTESIAN_FREEZE_WEIGHT"],
        },
    }

    @property
    def constraint_type(self) -> str:
        """hard coded constraint type"""
        return "backbone_freeze"

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
        
        csts[cidx] = ctor( atom_list, target_value, cst )

    return ResiduePairConstraint(
                                residue1,
                                residue2,
                                r1_atoms,
                                r2_atoms,
                                *csts
                                )

def create_backbone_freeze(stru: Structure) -> BackBoneFreeze:
    """constructor for BackboneFreeze. Constrain only C,CA,N"""
    atoms = stru.backbone_atoms()
    return BackBoneFreeze(atoms=atoms)

def create_distance_constraint() -> DistanceConstraint:
    """constructor for DistanceConstraint"""
# endregion
