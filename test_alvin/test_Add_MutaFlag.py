from Class_PDB import *

pdb_obj = PDB('./test_alvin/test_PDB/FAcD-FA-ASP.pdb') 

#def test_0():   
    #MutaFlag_0 = 'I202D'
    #pdb_obj.Add_MutaFlag(MutaFlag_0)
    #return pdb_obj.MutaFlags

def test_1():   
    MutaFlag_1 = 'I202D'
    pdb_obj.Add_MutaFlag(MutaFlag_1)
    assert pdb_obj.MutaFlags == [("I","A","202","D")]


def test_2():   
    MutaFlag_2 = ''
    pdb_obj.Add_MutaFlag(MutaFlag_2)
    assert pdb_obj.MutaFlags == [("")]