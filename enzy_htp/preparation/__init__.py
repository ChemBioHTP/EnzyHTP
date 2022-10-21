"""
preparation module for enzy_htp. Responsibilities include removal of waters, getting protonation state and generation of mutations.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""

from .pdb_line import PDBLine, read_pdb_lines
from .protonate import protonate_stru
from .clean import remove_solvent
