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
def test_pdbmd_without_job_manager_but_input_cpu_cores():
    test_dir_md = 'test/testfile_Class_PDB/MD_test/'
    # interface to PDBMD
    pdb_obj = PDB(f'{test_dir_md}FAcD_RA124M_ff.pdb', wk_dir=test_dir_md)
    pdb_obj.prmtop_path = f'{test_dir_md}FAcD_RA124M_ff.prmtop'
    pdb_obj.inpcrd_path = f'{test_dir_md}FAcD_RA124M_ff.inpcrd'
    # run MD
    with pytest.raises(TypeError) as e:
        pdb_obj.PDBMD(  engine='Amber_GPU', 
                        equi_cpu=1, 
                        if_cluster_job=0,
                        equi_cpu_cores=10
                         )
    assert 'ERROR: cpu_cores or equi_cpu_cores should be None if not submit to cluster.' in str(e.value)
    with pytest.raises(TypeError) as e:
        pdb_obj.PDBMD(  engine='Amber_CPU',
                        cpu_cores= 10,
                        if_cluster_job=0,
                        equi_cpu_cores=10
                         )
    assert 'ERROR: cpu_cores or equi_cpu_cores should be None if not submit to cluster.' in str(e.value)


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
    test_file_paths.extend(['mdcrd', 'mdinfo', 'submitted_job_ids.log'])
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
    test_file_paths.extend(['mdcrd', 'mdinfo', 'submitted_job_ids.log'])
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
    test_file_paths.extend(['mdcrd', 'mdinfo', 'submitted_job_ids.log'])
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
    test_file_paths.append('submitted_job_ids.log')
    for job in qm_jobs:
        test_file_paths.extend([job.sub_script_path, job.job_cluster_log])
    # assert result
    for out in qm_outs:
        assert os.path.isfile(out)
        assert os.path.getsize(out) != 0

@pytest.mark.qm
@pytest.mark.accre
def test_pdb2qmcluster_with_job_manager_str_res_setting():
    test_dir_qmcluster = 'test/testfile_Class_PDB/QMCluster_test/'
    # interface of PDB2QMCluster
    pdb_obj = PDB(f'{test_dir_qmcluster}FAcD_RA124M_ff.pdb', wk_dir=test_dir_qmcluster)
    pdb_obj.prmtop_path = f'{test_dir_qmcluster}FAcD_RA124M_ff.prmtop'
    pdb_obj.mdcrd = f'{test_dir_qmcluster}prod.mdcrd'
    pdb_obj.prepi_path = {'FAH':f'{test_dir_qmcluster}../ligands/ligand_FAH.prepin'}
    # arguments of PDB2QMCluster
    atom_mask = ':108,298'
    g_route = '# hf/6-31G(d) nosymm'
    res_str = '''#!/bin/bash
#SBATCH --partition=production
#SBATCH --job-name=EnzyHTP-test-str
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=30:00
#SBATCH --account=yang_lab_csb
'''
    with pytest.raises(TypeError) as e: # the case without nessessary information
        pdb_obj.PDB2QMCluster(
            atom_mask, 
            g_route=g_route, 
            if_cluster_job=1, 
            cluster=accre.Accre(), 
            job_array_size=20, 
            period=30,
            res_setting=res_str,
            cluster_debug=1
            )
    assert 'ERROR: Requires cpu_cores and cpu_mem input for configuring Gaussian input file if using CPU and provide res_setting not as a dict' in str(e.value)
    qm_outs = pdb_obj.PDB2QMCluster(
                atom_mask, 
                g_route=g_route, 
                if_cluster_job=1, 
                cluster=accre.Accre(), 
                job_array_size=20, 
                period=30,
                res_setting=res_str,
                cpu_cores=8,
                cpu_mem='2GB',
                cluster_debug=1
                )
    # track files
    qm_jobs = pdb_obj.qm_cluster_jobs
    qm_ins = list(map(lambda x: x.removesuffix('out')+'gjf',qm_outs))
    test_file_paths.extend(qm_outs)
    test_file_paths.extend(qm_ins)
    test_file_paths.append('submitted_job_ids.log')
    for job in qm_jobs:
        test_file_paths.extend([job.sub_script_path, job.job_cluster_log])
    # assert result
    for out in qm_outs:
        assert os.path.isfile(out)
        assert os.path.getsize(out) != 0

@pytest.mark.qm
def test_make_single_g16_job_loop():
    '''
    test the loop section in RunQM
    '''
    inp = [f'test_{i}.gjf' for i in range(10)]
    cluster = accre.Accre()
    jobs = []
    for gjf_path in inp:
        out_path = gjf_path.removesuffix('gjf')+'out'
        jobs.append(PDB._make_single_g16_job(gjf_path, out_path, cluster, PDB._get_default_res_setting_qmcluster(None)))
    for job in jobs:
        pass# assert len(job.sub_script_str.splitlines()) == len(jobs[0].sub_script_str.splitlines())

@pytest.mark.md
def test_pdbmin_local():
    """
    test for PDBMin running normally.
    need to apply an interactive shell to test
    """
    test_dir = 'test/testfile_Class_PDB/Min_test/'
    # interface to PDBMD
    pdb_obj = PDB(f'{test_dir}KE-07_ff.pdb', wk_dir=test_dir)
    pdb_obj.prmtop_path = f'{test_dir}KE-07.prmtop'
    pdb_obj.inpcrd_path = f'{test_dir}KE-07.inpcrd'
    # run MD
    pdb_obj.PDBMin(cycle=20000, engine='Amber_GPU', if_cluster_job=0)
    assert os.path.getsize(f'{test_dir}/cache/PDBMin/min.ncrst') != 0
    assert os.path.getsize(pdb_obj.path) != 0
    # clean up
    test_file_paths.extend([f'{test_dir}/cache/PDBMin/min.in',
                            f'{test_dir}/cache/PDBMin/min.out',
                            f'{test_dir}/cache/PDBMin/min.ncrst',
                            f'{test_dir}/KE-07_ff_min.pdb'])
    os.system(f"mv {pdb_obj.prmtop_path} {test_dir}KE-07.prmtop")
    os.system(f"mv {pdb_obj.inpcrd_path} {test_dir}KE-07.inpcrd")
    test_file_dirs.append(f'{test_dir}/cache/PDBMin')
    test_file_dirs.append(f'{test_dir}/cache')

@pytest.mark.md
def test_pdbmin_local_w_constrain():
    """
    test for PDBMin running with constrain applied. 
    need to apply an interactive shell to test
    """
    test_dir = 'test/testfile_Class_PDB/Min_test/'
    # interface to PDBMD
    pdb_obj = PDB(f'{test_dir}KE-07_ff.pdb', wk_dir=test_dir)
    pdb_obj.prmtop_path = f'{test_dir}KE-07.prmtop'
    pdb_obj.inpcrd_path = f'{test_dir}KE-07.inpcrd'
    # run MD
    # set up constrain
    pdb_obj.conf_min['DISANG']=pdb_obj.cache_path+'/PDBMin/0.rs'
    pdb_obj.conf_min['nmropt_rest'] = '1'
    pdb_obj.get_stru()
    a1=str(pdb_obj.stru.ligands[0].CAE.id)
    a2=str(pdb_obj.stru.ligands[0].H2.id)
    a3=str(pdb_obj.stru.chains[0].i101.OE2.id)
    pdb_obj.conf_min['rs_constraints']=[
                            {
                            'iat':[a2,a3],'r1': '0.50','r2': '0.96','r3': '2.50','r4': '3.50',
                            'rk2':'0.0','rk3':'100.0','ir6':'1','ialtd':'0'
                            },
                            {
                            'iat':[a1,a2,a3],'r1': '90.0','r2': '130.0','r3': '230.0','r4': '270.0',
                            'rk2':'100.0','rk3':'100.0','ir6':'1','ialtd':'0'
                            },
                        ]
    pdb_obj.PDBMin(cycle=20000, engine='Amber_GPU', if_cluster_job=0)
    assert os.path.getsize(f'{test_dir}/cache/PDBMin/min.ncrst') != 0
    assert os.path.getsize(pdb_obj.path) != 0
    # clean up
    test_file_paths.extend([f'{test_dir}/cache/PDBMin/min.in',
                            f'{test_dir}/cache/PDBMin/min.out',
                            f'{test_dir}/cache/PDBMin/min.ncrst',
                            f'{test_dir}/cache/PDBMin/0.rs'])
    test_file_paths.append(f'{test_dir}/KE-07_ff_min.pdb')
    os.system(f"mv {pdb_obj.prmtop_path} {test_dir}KE-07.prmtop")
    os.system(f"mv {pdb_obj.inpcrd_path} {test_dir}KE-07.inpcrd")
    test_file_dirs.append(f'{test_dir}/cache/PDBMin')
    test_file_dirs.append(f'{test_dir}/cache')

@pytest.mark.accre
def test_relax_with_rosetta():
    '''
    test relax_with_rosetta functions normally TODO
    '''
    test_dir = 'test/testfile_Class_PDB/ddg_test/'
    pdb_obj = PDB(f'{test_dir}KE-07.pdb', wk_dir=test_dir)
    pdb_obj.relax_with_rosetta(
        rosetta_home='/data/yang_lab/shaoqz/software/Rosetta313/main/',
        cluster=accre.Accre(),
        period=5)

def test_get_rosetta_lowest_score():
    '''
    test function works without abort TODO
    '''
    test_dir = 'test/testfile_Class_PDB/ddg_test/'
    target_idx = PDB.get_rosetta_lowest_score(f'{test_dir}score.sc')

@pytest.mark.accre
def test_get_rosetta_ddg():
    '''
    test get_rosetta_ddg works without abort
    '''
    test_dir = 'test/testfile_Class_PDB/ddg_test/'
    pdb_obj = PDB(f'{test_dir}KE-07.pdb', wk_dir=test_dir)
    ddg_results = pdb_obj.get_rosetta_ddg(
        rosetta_home='/data/yang_lab/Common_Software/Rosetta3.9/main/',
        muta_groups=[('D 7 I', 'L 2 I'), ('D 7 A')],
        relaxed_pdb=f'{test_dir}KE-07_relaxed.pdb',
        cluster=accre.Accre(),
        period=10)

def test_get_rosetta_ddg_result():
    '''
    test function works without abort
    '''
    test_dir = 'test/testfile_Class_PDB/ddg_test/'
    ddg = PDB.get_rosetta_ddg_result(f'{test_dir}example.ddg', 10)
    assert ddg.round(4) == 2.8953   

def test_get_cpptraj_rmsd_result():
    '''test function works without abort'''
    test_dir = 'test/testfile_Class_PDB/traj_ana_test/'
    result = PDB.get_cpptraj_rmsd_result(f'{test_dir}/rmsd.dat')
    assert result.round(4) == 0.3746

def test_get_sasa_ratio():
    '''test function works without abort'''
    test_dir = 'test/testfile_Class_PDB/traj_ana_test/'
    result = PDB.get_sasa_ratio(
        f'{test_dir}/FAcD_RA124M_ff.prmtop', 
        f'{test_dir}/prod.mdcrd',':1-297', 
        ':108,109,112,133,139,153,154,217,251,278', 
        ':298', 
        tmp_dir = test_dir)
    assert round(result, 4) == 6.2885

def test_get_ses_ratio():
    '''test function works without abort'''
    Config.debug = 1
    test_dir = 'test/testfile_Class_PDB/traj_ana_test/'
    result = PDB.get_ses_ratio(
        f'{test_dir}/FAcD_RA124M_ff.prmtop', 
        f'{test_dir}/prod.mdcrd',':1-297', 
        ':108,109,112,133,139,153,154,217,251,278', 
        ':298', 
        tmp_dir = test_dir)
    assert round(result, 4) == 0.0595

def test_get_residue_pka():
    '''test function works without abort. 
    Test getting pka for E101 of KE07R7'''
    test_pdb = 'test/testfile_Class_PDB/KE07R7.pdb'

    pdb_obj = PDB(test_pdb)
    result = pdb_obj.get_residue_pka(':101,222')
    assert result == [6.866247921750009, 9.868125271097515]

def test_make_mmpbsa_prmtops_use_ante_mmpbsa():
    '''test function works as expected'''
    ligand_mask = ':902'
    test_dir = 'test/testfile_Class_PDB/mmpbsa_test/'
    test_pdb = PDB(f"{test_dir}PuOrh_amber_aH_rmH_aH_ff.pdb", wk_dir=test_dir)
    test_pdb.prepi_path["FAD"] = f"{test_dir}../ligands/ligand_FAD.prepin"
    test_pdb.prepi_path["ACP"] = f"{test_dir}../ligands/ligand_ACP.prepin"
    test_pdb.frcmod_path["FAD"] = f"{test_dir}../ligands/ligand_FAD.frcmod"
    test_pdb.frcmod_path["ACP"] = f"{test_dir}../ligands/ligand_ACP.frcmod"
    test_pdb.prmtop = f'{test_dir}/PuOrh_amber_aH_rmH_aH.prmtop'
    Config.debug = 1

    assert test_pdb.make_mmpbsa_prmtops(ligand_mask) == (
        'test/testfile_Class_PDB/mmpbsa_test/temp/dr.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/dl.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/dc.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/sc.prmtop')

def test_make_mmpbsa_prmtops_no_ante_mmpbsa():
    '''test function works as expected'''
    ligand_mask = ':902'
    test_dir = 'test/testfile_Class_PDB/mmpbsa_test/'
    test_pdb = PDB(f"{test_dir}PuOrh_amber_aH_rmH_aH_ff.pdb", wk_dir=test_dir)
    test_pdb.prepi_path["FAD"] = f"{test_dir}../ligands/ligand_FAD.prepin"
    test_pdb.prepi_path["ACP"] = f"{test_dir}../ligands/ligand_ACP.prepin"
    test_pdb.frcmod_path["FAD"] = f"{test_dir}../ligands/ligand_FAD.frcmod"
    test_pdb.frcmod_path["ACP"] = f"{test_dir}../ligands/ligand_ACP.frcmod"
    test_pdb.prmtop = f'{test_dir}/PuOrh_amber_aH_rmH_aH.prmtop'
    Config.debug = 1

    assert test_pdb.make_mmpbsa_prmtops(ligand_mask, use_ante_mmpbsa=0) == (
        'test/testfile_Class_PDB/mmpbsa_test/temp/dr.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/dl.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/dc.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/sc.prmtop')
    
    test_file_paths.extend([
        'test/testfile_Class_PDB/mmpbsa_test/temp/dr.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/dl.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/dc.prmtop', 
        'test/testfile_Class_PDB/mmpbsa_test/temp/sc.prmtop'])
    test_file_dirs.append(f"{test_dir}temp/")


def test_run_mmpbsa():
    '''test function works as expected'''
    test_dir = 'test/testfile_Class_PDB/mmpbsa_test/'
    test_pdb = PDB(f"{test_dir}PuOrh_amber_aH_rmH_aH_ff.pdb", wk_dir=test_dir)
    dr_prmtop = f'{test_dir}data/dr.prmtop'
    dl_prmtop = f'{test_dir}data/dl.prmtop'
    dc_prmtop = f'{test_dir}data/dc.prmtop'
    sc_prmtop = f'{test_dir}data/sc.prmtop'
    traj_file = f'{test_dir}prod.mdcrd'
    answer_file = f'{test_dir}data/mmpbsa.dat'

    Config.debug = 1
    mmpbsa_out_file = test_pdb.run_mmpbsa(
        dr_prmtop, dl_prmtop, dc_prmtop, sc_prmtop, traj_file, 
        cluster=accre.Accre(),
        res_setting = {'account':'yang_lab'})
    
    assert os.path.getsize(mmpbsa_out_file) == os.path.getsize(answer_file)

    test_file_paths.append(mmpbsa_out_file)


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