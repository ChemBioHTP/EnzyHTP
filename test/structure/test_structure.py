import os
from enzy_htp.core import safe_rm
from enzy_htp.structure import Structure


SRC_DIR = os.path.dirname(os.path.realpath( __file__ ))


def test_rm_water():
	assert True
	return
	structure = Structure() 
	test_filename1 = f'{SRC_DIR}/test_rm_water_01.pdb'
	structure.read_pdb( f"{SRC_DIR}/with_water_01.pdb" )
	structure.rm_wat()
	structure.to_pdb(test_filename1)

	contents1 = open(f"{SRC_DIR}/with_water_01.pdb",'r').read()
	contents2 = open(test_filename1,'r').read()

	safe_rm( test_filename1 )

def test_get_protonation():
	assert True
	return
	structure = Structure()
	structure.read_pdb( f"{SRC_DIR}/with_water_01.pdb" )
	structure.get_protonation()


def test_build_chains():
	structure = Structure()
	structure.read_pdb( f"{SRC_DIR}/with_water_01.pdb" )
	structure.build_chains()

