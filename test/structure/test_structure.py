#TODO(CJ) Documentation
import os

from enzy_htp.structure import Structure, structure_from_pdb



CURR_DIR = os.path.dirname(os.path.abspath( __file__ ))
TEST_DIR = os.path.dirname( CURR_DIR )


def test_load_structure():
    """"""
    TEST_FILE = f"{TEST_DIR}/preparation/data/3NIR.pdb"
    structure_from_pdb(TEST_FILE)
    assert False
