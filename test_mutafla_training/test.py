import sys
#import pytest
#from absl import flags
#import run_enzyHTP
from Class_PDB import *
#from Bio.SeqUtils import IUPACData

pdb_path ='/scratch/ranx/EnzyHTP/Test_file/FAcD/FAcD-FA-ASP.pdb' 

pdb_file = PDB(pdb_path)

"""
def read_pdb_sequence(name):
    prev = '-1'
    prev_id = '-1'
    resi_id = []
    with open(name) as f :
        for line in f:
            toks = line.split()
            if len(toks)<1 : continue
            if toks[0] != 'ATOM': continue 
            if toks[3] != prev:
                seq.append(toks[3])
            prev = toks[3]
            if toks[4] != prev:
                resi_id.append(toks[4])
            prev_id = toks[4]
    return seq, resi_id


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

def test_4_fnx_add_mutaflag():
    mutaflag_4_test = 'ASN198Y'
    pdb_file.Add_MutaFlag(mutafla_3_test)
    assert pdb_file.MutaFlags == [("N", "A", "198","Y")]
