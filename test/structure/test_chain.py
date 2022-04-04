"""Testing the functionality of the enzy_htp.structure.Chain class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-20
"""
import os

from enzy_htp.structure import Chain, Residue


def test_proper_ctor_behavior():
    """Making sure that the default Chain() works."""
    chain = Chain('test', [])
    assert chain.name() == 'test'
    assert chain.empty()
    assert not chain.residues()

#TODO(CJ): add tests for residue accession and chain renaming
#TODO(CJ): need to add tests for checking if the chain has ceratin residues/residue types.
