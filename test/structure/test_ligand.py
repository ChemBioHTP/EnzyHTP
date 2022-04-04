"""Testing the Ligand() class in enzy_htp.structure. 

Author: Chris Jurich <cjurich2@huskers.unl.edu>
Date: 2022-04-03
"""
import pytest
import numpy as np

from enzy_htp.chemical import enum as renum
from enzy_htp.structure import Ligand


def test_constat_data():
    """Testing a variety of constant data methods that should work."""
    lig = Ligand("A.B.10", list())
    assert lig.is_ligand()
    assert not lig.is_canonical()
    assert not lig.is_metal()
    assert lig.rtype() == renum.ResidueType.LIGAND

def test_net_charge():
    """Checking that the net_charge attribute works properly for both default and set values."""
    lig = Ligand("A.B.10", list())
    assert lig.get_net_charge() is None
    lig = Ligand("A.B.10", list(), net_charge=1.0)
    assert np.isclose(lig.get_net_charge(), 1.0)


def test_build_bad_fname():
    """Ensuring the Ligand.build() method fails when a non .pdb files is supplied."""
    lig = Ligand("A.B.10", list())
    with pytest.raises(SystemExit) as exe:
         lig.build('test.txt')

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1 

