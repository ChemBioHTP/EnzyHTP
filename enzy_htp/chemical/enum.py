"""Enumerated types used to represent atomic and residue information in protein structures

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-02-20
"""
from enum import IntEnum


class ResidueType(IntEnum):
    CANONICAL = 0
    NONCANONICAL = 1
    UNKNOWN = 2
