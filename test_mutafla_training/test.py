import sys
#import pytest
#from absl import flags
#import run_enzyHTP
from Class_PDB import *

pdb_file = PDB('/scratch/ranx/EnzyHTP/Test_file/FAcD/FAcD-FA-ASP.pdb')

"""
flags.DEFINE_string('pdb_path_file',None, 'Path to the test pbd file.')
pdb_file = PDB('FLAGS.pdb_path_file')
FLAGS = flags.FLAGS
class TestClass:
    pdb_file = PDB('/scratch/ranx/EnzyHTP/Test_file/FAcD/FAcD-FA-ASP.pdb')

    def test1(self):
        self.pdb_file = 
"""


#@pytest.fixture(scope='function')
def test_1_fnx_add_mutaflag():
#    pdb_obj = PDB('/scratch/ranx/EnzyHTP/Test_file/FAcD/FAcD-FA-ASP.pdb')
    mutaflag_1_test = 'AA209C'
    pdb_file.Add_MutaFlag(mutaflag_1_test)
    assert pdb_file.MutaFlags == [("A","A", "209", "C")]

def test_2_fnx_add_mutaflag():
    mutaflag_2_test = 'A2557Y'
    pdb_file.Add_MutaFlag(mutaflag_2_test)
    assert pdb_file.MutaFlags == [("A","A","2557","C")]

def test_3_fnx_add_mutaflag():
    mutaflag_3_test = 'O198!'
    pdb_file.Add_MutaFlag(mutaflag_3_test)
    assert pdb_file.MutaFlags == [("O","A","198","H")]
