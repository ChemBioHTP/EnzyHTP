"""Testing the enumerated types found in enzy_htp.chemical.enum

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-04-03
"""
import pytest
from enzy_htp.chemical import enum as renum


def test_enum_value_comparisons():
    """Making sure that the enum values are correct."""
    assert renum.ResidueType.CANONICAL < renum.ResidueType.NONCANONICAL
    assert renum.ResidueType.NONCANONICAL < renum.ResidueType.SOLVENT
    assert renum.ResidueType.SOLVENT < renum.ResidueType.METAL
    assert renum.ResidueType.METAL < renum.ResidueType.LIGAND
    assert renum.ResidueType.LIGAND < renum.ResidueType.UNKNOWN

def test_converstion_to_integer():
    """Making sure conversion to integer works."""
    assert int(renum.ResidueType.CANONICAL) == 0
    assert int(renum.ResidueType.NONCANONICAL) == 1
    assert int(renum.ResidueType.SOLVENT) == 2
    assert int(renum.ResidueType.METAL) == 3
    assert int(renum.ResidueType.LIGAND) == 4
    assert int(renum.ResidueType.UNKNOWN) == 5

def test_conversion_from_integer():
    """Making sure conversion from integer works."""
    assert renum.ResidueType(0) == renum.ResidueType.CANONICAL
    assert renum.ResidueType(1) == renum.ResidueType.NONCANONICAL
    assert renum.ResidueType(2) == renum.ResidueType.SOLVENT
    assert renum.ResidueType(3) == renum.ResidueType.METAL
    assert renum.ResidueType(4) == renum.ResidueType.LIGAND
    assert renum.ResidueType(5) == renum.ResidueType.UNKNOWN


def test_residue_type_mapper():
    """Testing the RESIDUE_TYPE_MAPPER dict() found in enzy_htp.chemical.enum."""
    assert renum.RESIDUE_TYPE_MAPPER[renum.ResidueType(0)] == "CANONICAL"
    assert renum.RESIDUE_TYPE_MAPPER[renum.ResidueType(1)] == "NONCANONICAL"
    assert renum.RESIDUE_TYPE_MAPPER[renum.ResidueType(2)] == "SOLVENT"
    assert renum.RESIDUE_TYPE_MAPPER[renum.ResidueType(3)] == "METAL"
    assert renum.RESIDUE_TYPE_MAPPER[renum.ResidueType(4)] == "LIGAND"
    assert renum.RESIDUE_TYPE_MAPPER[renum.ResidueType(5)] == "UNKNOWN"

