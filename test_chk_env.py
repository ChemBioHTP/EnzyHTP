from helper import *
def test1():
    assert env_chk('pmemd') == 1
    assert env_chk('pmemd.MPI') == 1
    assert env_chk('g16') == 1
    assert env_chk('Multiwfn') == 1
def test2():
    assert env_chk('pmemd.cuda') == 0
def test3():
    assert env_chk(1) == -1
    assert env_chk([1]) == -1
    assert env_chk(['pmemd']) == -1
