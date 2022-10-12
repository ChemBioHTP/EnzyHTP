"""Testing the functionality implemented in enzy_htp.chemical.db.py. This file is a special case
as we do not want to update the database every time the unittests are run. As a result, we only 
test to see if the desired variables are defined within the database and that undefined variables
return 'None'.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-17
"""
import enzy_htp.chemical as chem


def test_variables_are_defined():
    """Testing that the required variables are defined."""
    assert chem.db.load_from_db("AA_LIST") is not None
    assert chem.db.load_from_db("IONIC_RADII") is not None
    assert chem.db.load_from_db("METAL_CENTER_MAP") is not None
    assert chem.db.load_from_db("METAL_MAPPER") is not None
    assert chem.db.load_from_db("ONE_LETTER_AA_MAPPER") is not None
    assert chem.db.load_from_db("RD_NON_LIGAND_LIST") is not None
    assert chem.db.load_from_db("RD_SOLVENT_LIST") is not None
    assert chem.db.load_from_db("RESIDUE_CATEGORIES") is not None
    assert chem.db.load_from_db("RESIDUE_CONNECTIVITY_MAP") is not None
    assert chem.db.load_from_db("RESIDUE_CONNECTIVITY_MAP_CTERMINAL") is not None
    assert chem.db.load_from_db("RESIDUE_CONNECTIVITY_MAP_NTERMINAL") is not None
    assert chem.db.load_from_db("RESIDUE_CTERMINAL_ATOM_LIST") is not None
    assert chem.db.load_from_db("RESIDUE_ELEMENT_MAP") is not None
    assert chem.db.load_from_db("RESIDUE_VOLUME_MAPPER") is not None
    assert chem.db.load_from_db("THREE_LETTER_AA_MAPPER") is not None
    assert chem.db.load_from_db("VDW_RADII") is not None
    assert chem.db.load_from_db("X_H_BOND_LENGTH") is not None


def test_undefined_variables_return_None():
    """Checking that undefined variables return the value None."""
    assert chem.db.load_from_db("dne") is None
    assert not chem.db.load_from_db("dne")
