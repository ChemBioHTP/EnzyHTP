# import sys

# sys.path.append('../EnzyHTP/')
from Class_PDB import *

def test_Add_MutaFlag_correct():
    # you need to create a pdb object first
    pdb_obj = PDB('./FAcD-FA-ASP.pdb') 
    MutaFlag_1 = 'I202D'

    pdb_obj.Add_MutaFlag(MutaFlag_1)

    assert pdb_obj.MutaFlags == [("I", "A", "202", "D")]

def test_Add_MutaFlag_wrong_type():
    # you need to create a pdb object first
    pdb_obj = PDB('./FAcD-FA-ASP.pdb') 
    MutaFlag_1 = 5

    pdb_obj.Add_MutaFlag(MutaFlag_1)

    assert pdb_obj.MutaFlags == []