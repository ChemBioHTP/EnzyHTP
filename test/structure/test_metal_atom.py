"""Testing the enzy_htp.structure.MetalUnit class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-04-03
"""
import os
import pytest
import numpy as np
from typing import List

import enzy_htp.chemical as chem
from enzy_htp.structure import (MetalUnit, PDBParser, Residue, Atom, residue_to_metal, PDBParser)

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"


def test_constant_data():
    """Making sure the constant data specialization for MetalUnit() is correct."""
    ma = MetalUnit(1, "ZN", [])

    assert not ma.is_ligand()
    assert ma.is_metal()
    assert ma.is_metal_center()
    assert not ma.is_canonical()
    assert ma.rtype == chem.ResidueType.METAL


def test_element():
    stru = PDBParser().get_structure(f"{DATA_DIR}1NVG.pdb")
    metalcenter: MetalUnit = stru.metalcenters[0]
    assert metalcenter.element == "Zn"


def test_get_donor_atoms():
    """test get donor atoms for a metal center from a stru"""
    stru = PDBParser().get_structure(f"{DATA_DIR}1NVG.pdb")
    metalcenter1: MetalUnit = stru.metalcenters[0]
    metalcenter2: MetalUnit = stru.metalcenters[1]

    assert len(metalcenter1.get_donor_atoms()) == 4
    assert len(metalcenter2.get_donor_atoms()) == 4


def test_get_donor_mapper():
    stru = PDBParser().get_structure(f"{DATA_DIR}1NVG.pdb")
    metalcenter1: MetalUnit = stru.metalcenters[0]
    donor_mapper = metalcenter1.get_donor_mapper()
    assert len(donor_mapper) == 4
    assert list(donor_mapper.values())[0][0].name == "OE2"


def test_get_radii_method_good_input():
    """Ensuring the MetalUnit.get_radii() method works correctly."""
    pdb_file = f"{DATA_DIR}/1NVG.pdb"
    metals: List[MetalUnit] = PDBParser().get_structure(pdb_file).metals
    assert np.isclose(metals[0].radius(), 0.88)
    assert np.isclose(metals[0].radius("vdw"), 1.39)

    assert np.isclose(metals[0].radius(), metals[1].radius())
    assert np.isclose(metals[0].radius("vdw"), metals[1].radius("vdw"))


def test_residue_to_metal():
    """Ensuring the residue_to_metal() method works."""
    residue = Residue(5, "ZN", list())
    metal: MetalUnit = residue_to_metal(residue)
    assert isinstance(metal, MetalUnit)
    assert metal.name == "ZN"


def test_clone():
    """Checking that the Ligand.clone() method returns a deepcopy of the current MetalUnit()."""
    pdb_file = f"{DATA_DIR}/1NVG.pdb"
    metals: List[MetalUnit] = PDBParser().get_structure(pdb_file).metals
    metal: MetalUnit = metals[0]
    metal_cpy: MetalUnit = metal.clone()
    assert isinstance(metal_cpy, MetalUnit)
    assert id(metal) != id(metal_cpy)
    for a1, a2 in zip(metal.atoms, metal_cpy.atoms):
        assert id(a1) != id(a2)
