from Class_PDB import *
obj = PDB('EnzyHTP/Test_file/FAcD/FAcD-FA-ASP.pdb')

def test():
    obj.Add_MutaFlag('A30T')
    assert obj.MutaFlags == [('A', 'A', '30', 'T')]
    obj.MutaFlags = []
def test2():
    obj.Add_MutaFlag('A11T')
    obj.Add_MutaFlag('A20T')
    assert obj.MutaFlags == [('A', 'A', '11', 'T'), ('A', 'A', '20', 'T')]
    obj.MutaFlags = []
def test3():
    obj.Add_MutaFlag('AA30R')
    assert obj.MutaFlags == [('A', 'A', '30', 'R')]
    obj.MutaFlags = []
def test4():
    obj.Add_MutaFlag('F31W')
    assert obj.MutaFlags[0] == ('F', 'A', '31', 'W')
    assert obj.MutaFlags[0][0] == 'F'
    assert obj.MutaFlags[0][1] == 'A'
    assert obj.MutaFlags[0][2] == '31'
    assert obj.MutaFlags[0][3] == 'W'
    obj.MutaFlags = []
def test5():
    obj.Add_MutaFlag(['V23T', 'W36F', 'V140C'])
    assert obj.MutaFlags == [('V', 'A', '23', 'T'), ('W', 'A', '36', 'F'), ('V', 'A', '140', 'C')]
