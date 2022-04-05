"""
structure module for enzy_htp. Describes hierarhical representation of Enzymes, including atoms, residues, chains and structures.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""

from .atom import Atom
from .ligand import Ligand
from .residue import Residue
from .solvent import Solvent, residue_to_solvent
from .metal_atom import MetalAtom, residue_to_metal
from .chain import Chain
from .structure import Structure, compare_structures, merge_right
from .structure_parser import structure_from_pdb, ligand_from_pdb
