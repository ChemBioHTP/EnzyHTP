from helper import *
import pytest

'''
Sans_check_test
'''
def test_pass():
    with pytest.raises(Exception):
        x = MutaFlag_san_check(['A','A',1,'C'])
        
def test_fail():
    with pytest.raises(Exception):
        x = MutaFlag_san_check(['A','A',1000000,'C'])
