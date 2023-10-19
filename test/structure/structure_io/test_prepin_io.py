"""Testing the PrepinParser class in the enzy_htp.structure.structure_io.prepin_io

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-10-17
"""
from pathlib import Path
from enzy_htp.structure.structure_io import PrepinParser

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/../data"

def test_prepin_parser_get_stru():
    """make sure function works as expected"""
    test_prepin = f"{DATA_DIR}/ligand_H5J.prepin"
    test_stru = PrepinParser().get_structure(test_prepin)
    print(test_stru)
