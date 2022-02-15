from Class_PDB import *

def run_Add_MutaFlag():
    # you need to create a pdb object first
    pdb_obj = PDB('./Test_file/FAcD/FAcD-FA-ASP.pdb')
    MutaFlag_1 = 'I202D'
    pdb_obj.Add_MutaFlag(MutaFlag_1)

    return pdb_obj, pdb_obj.MutaFlags, MutaFlag_1

def test_flag_is_string():
    #test if the given flag is a string
    a,b,c = run_Add_MutaFlag()
    assert type(c) == str

def test_false():
    #test if flag is only letters and numbers
    pdb_obj = PDB('./Test_file/FAcD/FAcD-FA-ASP.pdb')
    MutaFlag_1 = '?'
    assert MutaFlag_1.isalnum() == True
