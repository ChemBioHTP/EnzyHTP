"""Stores mappers and definitions for different types of solvents found in PDBS.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
"""
from typing import List

RD_NON_LIGAND_LIST: List[str] = ["CL", "EDO", "GOL"]
"""Common co-crystallized ligands in solvents."""

RD_SOLVENT_LIST: List[str] = ["HOH", "WAT"]
"""Common aliases for water."""