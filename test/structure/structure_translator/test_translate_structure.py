"""Testing correctness for the transate_structure() function found in enzy_htp.structure.structure_translator.api.py

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-05-27
"""
import pytest
import os

from enzy_htp.structure import (
    translate_structure,
    PDBParser,
    Structure
)

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"

def test_translate_structure_bad_naming_inputs() -> None:
    """Making sure translate_structure() errors for unsupported naming schemes."""
    val_error = None
    with pytest.raises(ValueError) as val_error:        
        translate_structure(None)
   
    assert val_error is not None

    val_error = None
    with pytest.raises(ValueError) as val_error:        
        translate_structure(None, start_naming="not_rosetta")
   
    assert val_error is not None
    
    val_error = None
    
    with pytest.raises(ValueError) as val_error:        
        translate_structure(None, end_naming="not_rosetta")

    assert val_error is not None

def test_translate_structure_correct_result() -> None:
    """Making sure the desired changes are made when using the translate_structure() function."""
    parser = PDBParser()
    stru:Structure = parser.get_structure(f"{DATA_DIR}/1Q4T_peptide_protonated.pdb")
    
    assert stru.residues[11].name == 'HID'
    
    translate_structure(stru, end_naming='rosetta')
    assert stru.residues[11].name == 'HIS'
    
    translate_structure(stru, start_naming='rosetta')
    assert stru.residues[11].name == 'HID'
