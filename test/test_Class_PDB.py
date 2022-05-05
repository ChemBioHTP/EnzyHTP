import copy
from glob import glob
from random import choice
import os
import pytest

from Class_PDB import PDB
from Class_Conf import Config
from core.clusters import accre
from helper import is_empty_dir
from AmberMaps import Resi_list, Resi_map2

test_file_paths = []
test_file_dirs = []
Config.debug=2
Config.Amber.conf_equi['nstlim'] = 50000
Config.Amber.conf_prod['nstlim'] = 500000

@pytest.mark.mutation
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
@pytest.mark.mutation
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

@pytest.mark.mutation
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

@pytest.mark.md
def test_PDB2FF_keep():
    pdb_obj = PDB('./test/testfile_Class_PDB/FAcD.pdb', wk_dir='./test/testfile_Class_PDB')
    prm_files = pdb_obj.PDB2FF(local_lig=1)
    test_file_paths.extend(prm_files) #clean up record
    test_file_paths.extend([pdb_obj.cache_path+'/leap.in', 
                            './tmp/tmp.inpcrd', 
                            pdb_obj.cache_path+'/leap.out',
                            pdb_obj.lig_dir+'/cache/ligand_temp2.pdb',
                            pdb_obj.lig_dir+'/cache/ligand_temp3.pdb',
                            pdb_obj.lig_dir+'/cache/ligand_temp.mol2',
                            './leap.log'])
    test_file_dirs.extend([pdb_obj.lig_dir+'/cache'])

    for f in prm_files:
        assert os.path.isfile(f)
        assert os.path.getsize(f) != 0 # prmtop will be 0K if failed

    assert len(pdb_obj.prepi_path) != 0

@pytest.mark.md
@pytest.mark.accre
def test_pdbmd_with_job_manager_capture_amber_err():
    test_dir_md = 'test/testfile_Class_PDB/MD_test/'
    # interface to PDBMD
    pdb_obj = PDB(f'{test_dir_md}FAcD_RA124M_ff.pdb', wk_dir=test_dir_md)
    pdb_obj.prmtop_path = f'{test_dir_md}FAcD_RA124M_ff.prmtop'
    pdb_obj.inpcrd_path = f'{test_dir_md}FAcD_RA124M_ff.inpcrd'
    # run MD
    with pytest.raises(Exception) as e:
        pdb_obj.PDBMD(  engine='Amber_GPU', 
                        equi_cpu=0, 
                        if_cluster_job=1,
                        cluster=accre.Accre(),
                        period=30,
                        res_setting={'account':'csb_gpu_acc'},
                        cluster_debug=1 )
        assert 'STOP PMEMD Terminated Abnormally!' in str(e.value)
    # track files
    test_file_paths.extend(['mdcrd', 'mdinfo'])
    test_file_paths.extend(glob(f'{test_dir_md}MD/*'))
    test_file_dirs.extend([f'{test_dir_md}cache',f'{test_dir_md}MD'])

@pytest.mark.md
@pytest.mark.accre
def test_pdbmd_with_job_manager_no_equi_cpu():
    test_dir_md = 'test/testfile_Class_PDB/MD_test_full_GPU/'
    # interface to PDBMD
    pdb_obj = PDB(f'{test_dir_md}GPU_test_ff.pdb', wk_dir=test_dir_md)
    pdb_obj.prmtop_path = f'{test_dir_md}GPU_test_ff.prmtop'
    pdb_obj.inpcrd_path = f'{test_dir_md}GPU_test_ff.inpcrd'
    # run MD
    nc_path = pdb_obj.PDBMD(engine='Amber_GPU', 
                            equi_cpu=0, 
                            if_cluster_job=1,
                            cluster=accre.Accre(),
                            period=30,
                            res_setting={'account':'csb_gpu_acc'},
                            cluster_debug=1 )
    assert os.path.isfile(nc_path)
    assert os.path.getsize(nc_path) != 0
    # track files
    os.remove(nc_path)
    test_file_paths.extend(['mdcrd', 'mdinfo'])
    test_file_paths.extend(glob(f'{test_dir_md}MD/*'))
    test_file_dirs.extend([f'{test_dir_md}cache',f'{test_dir_md}MD'])

@pytest.mark.md
@pytest.mark.accre
def test_pdbmd_with_job_manager_equi_cpu():
    test_dir_md = 'test/testfile_Class_PDB/MD_test/'
    # interface to PDBMD
    pdb_obj = PDB(f'{test_dir_md}FAcD_RA124M_ff.pdb', wk_dir=test_dir_md)
    pdb_obj.prmtop_path = f'{test_dir_md}FAcD_RA124M_ff.prmtop'
    pdb_obj.inpcrd_path = f'{test_dir_md}FAcD_RA124M_ff.inpcrd'
    # run MD
    nc_path = pdb_obj.PDBMD(engine='Amber_GPU', 
                            equi_cpu=1, 
                            if_cluster_job=1,
                            cluster=accre.Accre(),
                            period=30,
                            res_setting={'account':'csb_gpu_acc'},
                            res_setting_equi_cpu={'account':'yang_lab_csb'},
                            cluster_debug=1 )
    assert os.path.isfile(nc_path)
    assert os.path.getsize(nc_path) != 0
    os.remove(nc_path)
    # track files
    test_file_paths.extend(['mdcrd', 'mdinfo'])
    test_file_paths.extend(glob(f'{test_dir_md}MD/*'))
    test_file_dirs.extend([f'{test_dir_md}cache',f'{test_dir_md}MD'])

@pytest.mark.qm
def test_get_default_res_setting_qmcluster():
    config_default_before = copy.deepcopy(Config.Gaussian.QMCLUSTER_CPU_RES)
    correct_res_setting = copy.deepcopy(config_default_before)
    correct_res_setting['node_cores'] = '8'
    assert PDB._get_default_res_setting_qmcluster({'node_cores': '8'}) == correct_res_setting
    assert PDB._get_default_res_setting_qmcluster('''test_str''') == '''test_str'''
    assert Config.Gaussian.QMCLUSTER_CPU_RES == config_default_before # should not change the global default

@pytest.mark.qm
@pytest.mark.accre
def test_pdb2qmcluster_with_job_manager():
    test_dir_qmcluster = 'test/testfile_Class_PDB/QMCluster_test/'
    # interface of PDB2QMCluster
    pdb_obj = PDB(f'{test_dir_qmcluster}FAcD_RA124M_ff.pdb', wk_dir=test_dir_qmcluster)
    pdb_obj.prmtop_path = f'{test_dir_qmcluster}FAcD_RA124M_ff.prmtop'
    pdb_obj.mdcrd = f'{test_dir_qmcluster}prod.mdcrd'
    pdb_obj.prepi_path = {'FAH':f'{test_dir_qmcluster}../ligands/ligand_FAH.prepin'}
    # arguments of PDB2QMCluster
    atom_mask = ':108,298'
    g_route = '# hf/6-31G(d) nosymm'
    qm_outs = pdb_obj.PDB2QMCluster(
                        atom_mask, 
                        g_route=g_route, 
                        if_cluster_job=1, 
                        cluster=accre.Accre(), 
                        job_array_size=20, 
                        period=30,
                        res_setting={'account':'yang_lab_csb'},
                        cluster_debug=1
                        )
    # track files
    qm_jobs = pdb_obj.qm_cluster_jobs
    qm_ins = list(map(lambda x: x.removesuffix('out')+'gjf',qm_outs))
    test_file_paths.extend(qm_outs)
    test_file_paths.extend(qm_ins)
    for job in qm_jobs:
        test_file_paths.extend([job.sub_script_path, job.job_cluster_log])
    # assert result
    for out in qm_outs:
        assert os.path.isfile(out)
        assert os.path.getsize(out) != 0

### utilities ###
@pytest.mark.clean
def test_clean_files():
    # clean files
    for f_path in test_file_paths:
        if os.path.isfile(f_path):
            os.remove(f_path)
            print(f'removed {f_path}')
        else:
            print(f'no such file {f_path}')
            
@pytest.mark.clean            
def test_clean_dirs():
    # clean dirs
    for f_dir in test_file_dirs:
        if is_empty_dir(f_dir) != 2:
            if is_empty_dir(f_dir):
                os.removedirs(f_dir)
                print(f'removed {f_dir}')
            else:
                print(f'{f_dir} is not empty. wont rm it')