from Class_PDB import *

def test_Add_MutaFlag():
    pdb_obj = PDB('./test/testfile_Class_PDB/FAcD.pdb') # you need to create a pdb object first
    MutaFlag_1 = 'I202D'
    pdb_obj.Add_MutaFlag(MutaFlag_1)

    assert pdb_obj.MutaFlags == [('I', 'A', '202', 'D')]