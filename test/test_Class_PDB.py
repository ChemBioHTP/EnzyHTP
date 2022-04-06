from random import choice
from Class_PDB import PDB
from AmberMaps import Resi_map2


def _random_gen_good_MutaFlag_for_test(stru, abbr=0):
    '''
    this should be seperate from Add_MutaFlag
    '''
    if abbr:
        chain = stru.chains[0]
    else:
        chain = choice(stru.chains)
    residue = choice(chain.residues)
    chain_id = chain.id
    residue_name = Resi_map2[residue.name]
    residue_id = str(residue.id)
    target_residue_name = choice(Resi_list[:-1])

    if abbr:
        test_flag = residue_name + residue_id + target_residue_name
    else:
        test_flag = residue_name + chain_id + residue_id + target_residue_name
    
    correct_answer = (residue_name, chain_id, residue_id, target_residue_name)

    return test_flag, correct_answer


# how should we determine how many types of input should one function contain.
def test_Add_MutaFlag_good_assign_mutaflag_canonical():
    '''
    test assign a specific reasonable mutation to a canonical amino acid
    '''
    pdb_obj = PDB('./test/testfile_Class_PDB/FAcD.pdb')  # what should we do if there's a general input for all test?
    # create an random input of "XA##B"
    pdb_obj.get_stru()

    # test single input
    test_flag_0, correct_answer = _random_gen_good_MutaFlag_for_test(pdb_obj.stru)
    print('test_Add_MutaFlag_assign_good_mutaflag_canonical(): MutaFlag used for test: ' + test_flag_0)
    pdb_obj.Add_MutaFlag(test_flag_0)
    assert pdb_obj.MutaFlags == [correct_answer]

    pdb_obj.MutaFlags = []
    # test list input
    test_flag_list = []
    correct_answer = []
    for i in range(10):  # pylint: disable:=unused-variable
        flag_i, answer_i = _random_gen_good_MutaFlag_for_test(pdb_obj.stru)
        test_flag_list.append(flag_i)
        correct_answer.append(answer_i)
    print('test_Add_MutaFlag_assign_good_mutaflag_canonical(): MutaFlag used for test:', test_flag_list)
    pdb_obj.Add_MutaFlag(test_flag_list)
    assert pdb_obj.MutaFlags == correct_answer

    pdb_obj.MutaFlags = []
    # test abbreviated input
    test_flag_1, correct_answer = _random_gen_good_MutaFlag_for_test(pdb_obj.stru, abbr=1)
    print('test_Add_MutaFlag_assign_good_mutaflag_canonical(): MutaFlag used for test: ' + test_flag_1)
    pdb_obj.Add_MutaFlag(test_flag_1)
    assert pdb_obj.MutaFlags == [correct_answer]

def test_Add_MutaFlag_good_random():
    '''
    test random generate a mutation
    '''
    pdb_obj = PDB('./test/testfile_Class_PDB/FAcD.pdb')  # what should we do if there's a general input for all test?

    pdb_obj.Add_MutaFlag('r')

    assert pdb_obj.MutaFlags # How should we test this? It's just write the same thing

# good non-canonical
# bad input wrong original residue
# bad input chain index out of range
# bad input residue index out of range
# bad input self mutation
# very-bad input random character
# very-bad input other obj type

# good PDB multichains
# good PDB ligand
# good PDB non-c AA
# good PDB solvent
# bad PDB

def test_select_residues_by_distance():
    pdb_obj = PDB('./test/testfile_Class_PDB/FAcD.pdb')
    
    #return pdb_obj.select_residues_geo('dis', ['cnt',-6.144,-65.450,84.041], [10.0,15.0], ['PRO', 'SER'])
    #return pdb_obj.select_residues_geo('dis', ['pro',0,3,'CA'], [10.0,15.0], ['PRO', 'SER'])
    #return pdb_obj.select_residues_geo('dis', ['lig',0,'F'], [10.0,15.0], ['PRO', 'SER'])
    #return pdb_obj.select_residues_geo('dis', ['met',3], [10.0,15.0], ['PRO', 'SER'])

    return pdb_obj.select_residues_geo('dis', ['met',3], [10.0,15.0], ['all'])
print(test_select_residues_by_distance())
