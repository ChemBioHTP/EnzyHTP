"""Testing the enzy_htp.structure.MetalAtom class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-04-03
"""
import os
import pytest
import numpy as np
from biopandas.pdb import PandasPdb

import enzy_htp.chemical as chem
from enzy_htp.structure import MetalAtom, structure_from_pdb


CURR_DIR = os.path.dirname(os.path.abspath( __file__ ))
DATA_DIR = f"{CURR_DIR}/data/"

def test_constant_data():
    """Making sure the constant data specialization for MetalAtom() is correct."""
    ma = MetalAtom('A.A.1', list())

    assert not ma.is_ligand()
    assert ma.is_metal()
    assert ma.is_metal_center()
    assert not ma.is_canonical()
    assert ma.rtype() == chem.ResidueType.METAL 


def test_get_radii_method_good_input():
    """Ensuring the MetalAtom.get_radii() method works correctly."""
    pdb_file = f"{DATA_DIR}/1NVG.pdb"
    metals = list(filter(lambda r: r.is_metal(), structure_from_pdb( pdb_file ).chains[0].residues()))
    print(metals[0].atom_name_)
    assert np.isclose(metals[0].get_radii(),0.88)
    assert np.isclose(metals[0].get_radii('vdw'),1.39)

    assert np.isclose(metals[0].get_radii(),metals[1].get_radii())
    assert np.isclose(metals[0].get_radii('vdw'),metals[1].get_radii('vdw'))


def test_residue_to_metal():
    """Ensuring the residue_to_metal() method works."""
    assert False
