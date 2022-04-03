from helper import MutaFlag_decode as decode
import pytest

def test_regular_input():
    x = decode('XA11Y')
    assert x == ('X', 'A', '11', 'Y')
    assert x[0].isalpha()
    assert x[1].isalpha()
    assert x[2].isdigit()
    assert x[3].isalpha()

def test_without_ChainIndex():
    x = decode('X11Y')
    assert x == ('X', 'A', '11', 'Y')
    assert x[0].isalpha()
    assert x[1].isalpha()
    assert x[2].isdigit()
    assert x[3].isalpha()

def test_unregulated_format():
    x = decode('  X a11     1y ')
    assert x == ('X', 'A', '111', 'Y')
    assert x[0].isalpha()
    assert x[1].isalpha()
    assert x[2].isdigit()
    assert x[3].isalpha()

def test_error_First3Digit():
    with pytest.raises(TypeError): 
        assert decode('  X 1a1     1A ')

def test_error_WrongInput():
    with pytest.raises(TypeError): 
        assert decode('  X 1     1A1 ')
