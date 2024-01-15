"""Testing the Mol2Parser class in the enzy_htp.structure.structure_io.mol2_io.py

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-10-23
"""
from pathlib import Path

import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_io import Mol2Parser

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/../data"

def test_mol2_parser_get_stru():
    """make sure function works as expected"""
    test_prepin = f"{DATA_DIR}/TYQ.mol2"
    test_stru = Mol2Parser().get_structure(test_prepin)
    assert len(test_stru.residues) == 1
