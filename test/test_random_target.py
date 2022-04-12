from random import choice
from Class_PDB import *
from AmberMaps import Resi_map2


def test_random_mut_generator():
    facd_test_obj = Class_PDB.PDB('../Test_file/FAcD/FAcD_FA_ASP.pdb')
    facd_test_obj.random_mut_generator(12, 'A')


print(test_random_mut_generator())
