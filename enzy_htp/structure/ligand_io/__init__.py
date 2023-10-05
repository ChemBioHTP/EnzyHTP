"""The I/O submodule for generation/construction of Ligand objects from files and other data structures. The Ligand() class uniquely needs 
I/O methods as bond information is frequently needed for these molecules given their customized topologies.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-10-05
"""

from .mol2_io import Mol2Parser
