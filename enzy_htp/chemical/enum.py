"""Enumerated types used to represent atomic and residue information in protein structures

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-02-20
"""
from enum import IntEnum


class ResidueType(IntEnum):
    """Enumerated type representing types of residues found in PDB entries."""
    CANONICAL = 0
    NONCANONICAL = 1
    SOLVENT = 2
    METAL = 3
    LIGAND = 4
    UNKNOWN = 5


RESIDUE_TYPE_MAPPER = {
    ResidueType.CANONICAL: "CANONICAL",
    ResidueType.NONCANONICAL: "NONCANONICAL",
    ResidueType.SOLVENT: "SOLVENT",
    ResidueType.METAL: "METAL",
    ResidueType.LIGAND: "LIGAND",
    ResidueType.UNKNOWN: "UNKNOWN"
}
