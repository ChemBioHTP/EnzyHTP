"""
structure module for enzy_htp. Describes hierarhical representation of Enzymes, including atoms, residues, chains and structures.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""

from .atom import Atom
from .residue import Residue
from .structure import Structure, structure_from_pdb
from .pdb_line import PDBLine, read_pdb_lines
from .pdb_prepper import PDBPrepper
