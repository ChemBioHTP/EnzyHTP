"""Testing the PrepinParser class in the enzy_htp.structure.structure_io.prepin_io

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-10-17
"""
from pathlib import Path

import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_io import PrepinParser

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/../data"

def test_prepin_parser_get_stru():
    """make sure function works as expected"""
    test_prepin = f"{DATA_DIR}/ligand_H5J.prepin"
    test_stru = PrepinParser().get_structure(test_prepin)
    print(test_stru)

def test_deduce_coord_end():
    """test using an example file"""
    test_prepin = f"{DATA_DIR}/ligand_H5J.prepin"
    assert PrepinParser._deduce_coord_end(fs.lines_from_file(test_prepin)) == 26

def test_parse_iopr():
    """test using an example file"""
    test_prepin = f"{DATA_DIR}/ligand_H5J.prepin"
    assert PrepinParser._parse_iopr(fs.lines_from_file(test_prepin)[26:]) == {
        "LOOP": [["CAC", "CAI"],
                 ["CAK", "CAJ"]],
        "IMPROPER": [
            ["CAI", "OAB", "NAL", "OAA"],
            ["CAC", "CAF", "CAI", "NAL"],
            ["CAI", "CAJ", "CAF", "H1"],
            ["CAK", "CAF", "CAJ", "CAE"],
            ["CAJ", "H2", "CAE", "NAG"],
            ["CAD", "CAJ", "CAK", "OAH"],
            ["CAK", "CAC", "CAD", "H3"],
            ["CAD", "CAI", "CAC", "H4"]]}

def test_parse_prepin_file():
    """test using an example file. Simply make sure no error occurs."""
    test_prepin = f"{DATA_DIR}/ligand_H5J.prepin"
    PrepinParser._parse_prepin_file(test_prepin)
