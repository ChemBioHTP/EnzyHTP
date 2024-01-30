"""
structure module for enzy_htp. Describes hierarhical representation of Enzymes, including atoms, 
residues, chains and structures in a doubly-linked manner

Structural description utilizes polymorphic specialization of base Residue() class into Ligand(),
MetalUnit() and Solvent() derived classes.

Structures are loaded through the PDBParser().get_structure() method defined in enzy_htp.structure_io.pdb_io.py

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: QZ Shao <shaoqz@icloud.com>
Date: 2022-03-19
"""

from .atom import Atom
from .ligand import Ligand
from .residue import Residue
from .solvent import Solvent, residue_to_solvent
from .metal_atom import MetalUnit, residue_to_metal
from .modified_residue import ModifiedResidue, residue_to_modified_residue
from .noncanonical_base import NonCanonicalBase
from .chain import Chain
from .structure import Structure
from .structure_ensemble import StructureEnsemble
from .structure_io import PDBParser, Mol2Parser, PrepinParser

from .structure_region import (
    StructureRegion,
    create_region_from_selection_pattern,
    create_region_from_residue_keys,
    create_region_from_full_stru,
)

from .structure_constraint import (
    StructureConstraint,
    CartesianFreeze,
    DistanceConstraint,
    AngleConstraint,
    DihedralConstraint,
    ResiduePairConstraint,
    structure_constraints_from_xml,
    create_residue_pair_constraint,
    create_cartesian_freeze,
    create_backbone_freeze,
    create_distance_constraint,
    create_angle_constraint,
    create_dihedral_constraint,
    merge_cartesian_freeze,
    freeze_hydrogen_bonds 
)

from .structure_selection_class import StruSelection
