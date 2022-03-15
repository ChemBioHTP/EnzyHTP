from Class_PDB import *
from helper import MutaFlag_decode
import pytest


def test_str():
    with pytest.raises(Exception):
        MutaFlag_decode(['XA11Y'])
        assert resi_1 == 'X'
        assert chain_id == 'A'
        assert resi_id == '11'
        assert resi_2 == 'Y'
    
def test_str_no_chainid():
    with pytest.raises(Exception):
        MutaFlag_decode(['X11Y'])
        assert resi_1 == 'X'
        assert chain_id == 'A'
        assert resi_id == '11'
        assert resi_2 == 'Y'