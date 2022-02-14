import os
from enzy_htp.core import safe_rm
from enzy_htp.structure import Structure


SRC_DIR = os.path.dirname(os.path.realpath( __file__ ))

def pdbs_identical( pdb1, pdb2 ):
	pass

def test_rm_water():
	structure = Structure()
	test_filename1 = f'{SRC_DIR}/test_rm_water_01.pdb'
	structure.read_pdb( f"{SRC_DIR}/with_water_01.pdb" )
	structure.rm_wat()
	structure.to_pdb(test_filename1)

	contents1 = open(f"{SRC_DIR}/with_water_01.pdb",'r').read()
	contents2 = open(test_filename1,'r').read()

	safe_rm( test_filename1 )
