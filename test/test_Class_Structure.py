import pytest

from Class_Structure import (
    Structure,
    Chain,
    Residue,
    Atom,
)

DATA_DIR = 'test/testfile_Class_Structure/'

def test_structure_get_connect():
    '''test if the function act as expected'''
    stru = Structure.fromPDB(f'{DATA_DIR}maa_test/3cfr-slp-pea_aH.pdb')
    stru.get_connect(prepi_path={
        'LLP' : f'{DATA_DIR}maa_test/llp.prepin',
        'SLP' : f'{DATA_DIR}maa_test/ligand_SLP.prepin'
    })
    # print(*[(str(atom), atom.connect) for atom in stru.find_name_residue('LLP')[0]], sep='\n')
    llp = stru.find_name_residue('LLP')[0]
    assert len(llp.N.connect) == 3
    assert len(llp.CA.connect) == 4
    assert len(llp.C.connect) == 3
    slp = stru.find_name_residue('SLP')[0]
    assert len(slp.C37.connect) == 3
