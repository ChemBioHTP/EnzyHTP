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
    assert test_stru.residues[0]

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

def test_build_atoms():
    """test using an example file."""
    test_prepin = f"{DATA_DIR}/ligand_H5J.prepin"
    prepin_data = PrepinParser._parse_prepin_file(test_prepin)
    # convert to cartesian coordinate
    PrepinParser._deduce_cartesian_coord(prepin_data)

    # build Structure()
    atoms = PrepinParser._build_atoms(prepin_data)
    assert len(atoms) == 16

def test_connect_atoms():
    """test using an example file."""
    test_prepin = f"{DATA_DIR}/ligand_H5J.prepin"
    prepin_data = PrepinParser._parse_prepin_file(test_prepin)
    # convert to cartesian coordinate
    PrepinParser._deduce_cartesian_coord(prepin_data)

    # build Structure()
    atoms = PrepinParser._build_atoms(prepin_data)
    PrepinParser._connect_atoms(atoms, prepin_data)

    answer = [
        ("OAB",  ["NAL"]),
        ("NAL",  ["OAB", "OAA", "CAI"]),
        ("OAA",  ["NAL"]),
        ("CAI",  ["NAL", "CAF", "CAC"]),
        ("CAF",  ["CAI", "H1", "CAJ"]),
        ("H1",   ["CAF"]),
        ("CAJ",  ["CAF", "CAE", "CAK"]),
        ("CAE",  ["CAJ", "H2", "NAG"]),
        ("H2",   ["CAE"]),
        ("NAG",  ["CAE", "OAH"]),
        ("OAH",  ["NAG", "CAK"]),
        ("CAK",  ["OAH", "CAD", "CAJ"]),
        ("CAD",  ["CAK", "H3", "CAC"]),
        ("H3",   ["CAD"]),
        ("CAC",  ["CAD", "H4", "CAI"]),
        ("H4",   ["CAC"]),]

    for atom, answer_atom in zip(atoms, answer):
        assert (atom.name, [catom.name for catom in atom.connect_atoms]) == answer_atom

def test_connect_atoms_noncyclic():
    """test using an example file."""
    test_prepin = f"{DATA_DIR}/PUT_any.prepin"
    prepin_data = PrepinParser._parse_prepin_file(test_prepin)
    # convert to cartesian coordinate
    PrepinParser._deduce_cartesian_coord(prepin_data)

    # build Structure()
    atoms = PrepinParser._build_atoms(prepin_data)
    PrepinParser._connect_atoms(atoms, prepin_data)

    answer = [
        ("N5", ['H15', 'H16', 'C2']),
        ("H15",  ['N5']),
        ("H16",  ['N5']),
        ("C2", ['N5', 'H8', 'H9', 'C1']),
        ("H8", ['C2']),
        ("H9", ['C2']),
        ("C1", ['C2', 'H7', 'H12', 'C3']),
        ("H7", ['C1']),
        ("H12",  ['C1']),
        ("C3", ['C1', 'H10', 'H11', 'C4']),
        ("H10",  ['C3']),
        ("H11",  ['C3']),
        ("C4", ['C3', 'H13', 'H14', 'N6']),
        ("H13",  ['C4']),
        ("H14",  ['C4']),
        ("N6", ['C4', 'H17', 'H18', 'H19']),
        ("H17",  ['N6']),
        ("H18",  ['N6']),
        ("H19",  ['N6']),]

    for atom, answer_atom in zip(atoms, answer):
        assert (atom.name, [catom.name for catom in atom.connect_atoms]) == answer_atom
