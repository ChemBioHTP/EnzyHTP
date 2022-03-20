""" Testing for the enzy_htp.preparation.PDBLine() class

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os
import pytest

from enzy_htp.preparation import PDBLine, read_pdb_lines
from enzy_htp.core import UnsupportedFileType, lines_from_file


TEST_DIR = os.path.dirname(os.path.abspath( __file__ ))


def test_basic_getters():
    """Ensuring that the basic getters for PDB line record type function."""

    end_line="END                                                                             "
    end_line=PDBLine(end_line)
    assert end_line.is_END()

    ter_line="TER     893      ASN A  46                                                      "
    ter_line=PDBLine(ter_line)
    assert ter_line.is_TER()

    hetatm_line="HETATM  894  C1  EOH A2001      18.245  -4.085  12.277  0.60 10.06           C  "
    hetatm_line=PDBLine(hetatm_line)
    assert hetatm_line.is_HETATM()

    atom_line="ATOM    892 HD22 ASN A  46      15.661  -2.840  13.786  1.00  7.87           H  "
    atom_line=PDBLine(atom_line)
    assert atom_line.is_ATOM()

def test_read_pdb_lines():
    """Checking that read_pdb_lines() works for both good and bad input files."""
    
    with pytest.raises(UnsupportedFileType) as exe:
        read_pdb_lines("not-a-pdb.txt")

    assert exe
    assert exe.type == UnsupportedFileType

    pdb_file = f"{TEST_DIR}/data/3NIR.pdb"
    assert read_pdb_lines( pdb_file )


def test_round_trip():
    """Testing that the round trip of a PDB lines actually works."""
    pdb_file = f"{TEST_DIR}/data/3NIR.pdb"
    raw_pdb_lines = lines_from_file( pdb_file )
    pdb_lines =  read_pdb_lines( pdb_file )
    rt_pdb_lines = list(map(lambda pl : pl.build(), pdb_lines ))
    for a,b in zip(raw_pdb_lines, rt_pdb_lines):
        print(f"'{a}'")
        print(f"'{b}'")
        assert a == b

