"""
structure module for enzy_htp. Describes hierarhical representation of Enzymes, including atoms, 
residues, chains and structures in a doubly-linked manner

Structural description utilizes polymorphic specialization of base Residue() class into Ligand(),
MetalUnit() and Solvent() derived classes.

Structures are loaded through the PDBParser().get_structure() method defined in enzy_htp.structure_io.pdb_io.py

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""

from .atom import Atom
from .ligand import Ligand
from .residue import Residue
from .solvent import Solvent, residue_to_solvent
from .metal_atom import MetalUnit, residue_to_metal
from .chain import Chain
from .structure import Structure
from .structure_io import PDBParser
from .structure_parser import structure_from_pdb, ligand_from_pdb  # TODO remove later
