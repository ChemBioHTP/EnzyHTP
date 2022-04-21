import pytest

from core import clusters
from core.job_manager import *

command_2_run = ['g16 < xxx.gjf > xxx.out']
env_settings_list =  [  'module load GCC/6.4.0-2.28  OpenMPI/2.1.1', 
                        'module load Amber/17-Python-2.7.14', 
                        'module unload Python/2.7.14 numpy/1.13.1-Python-2.7.14']
res_keywords_dict = {   'core_type' : 'cpu',
                        'node_cores' : '24',
                        'job_name' : 'job_name',
                        'partition' : 'production',
                        'mem_per_core' : '4G',
                        'walltime' : '24:00:00',
                        'account' : 'xxx'}
env_settings_str = '''export AMBERHOME=/dors/csb/apps/amber19/
export CUDA_HOME=$AMBERHOME/cuda/10.0.130
export LD_LIBRARY_PATH=$AMBERHOME/cuda/10.0.130/lib64:$AMBERHOME/cuda/RHEL7/10.0.130/lib:$LD_LIBRARY_PATH'''
res_keywords_str = '''#!/bin/bash
#SBATCH --job-name=ETIY127A          # Assign an 8-character name to your job
#SBATCH --account=xxx
#SBATCH --partition=turing
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=6
#SBATCH --mem=50G
#SBATCH --time=5-00:00:00         # Total run time limit (HH:MM:SS)
#SBATCH --export=ALL
'''
sub_script_str = '''#!/bin/bash
#SBATCH --partition=debug
#SBATCH --job-name=JM-test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=30:00
#SBATCH --account=yang_lab_csb

module load Gaussian/16.B.01
mkdir $TMPDIR/$SLURM_JOB_ID
export GAUSS_SCRDIR=$TMPDIR/$SLURM_JOB_ID

g16 < TS-2-dp-opt.gjf > TS-2-dp-opt.out'''

cluster = accre.Accre()

def test_ClusterJob_submit_by_list():
    job = ClusterJob.config_job(
        commands = command_2_run,
        cluster = cluster,
        env_settings = env_settings_list,
        res_keywords = res_keywords_dict
    )
    job.submit(debug=1, script_path='./test/core/test_file/test.cmd')
    #assert job.submit(debug=1, script_path='./test/core/test_file/test.cmd') == 'sbatch ./test/core/test_file/test.cmd'

def test_ClusterJob_submit_by_str():
    job = ClusterJob.config_job(
        commands = command_2_run,
        cluster = cluster,
        env_settings = env_settings_str,
        res_keywords = res_keywords_str
    )
    job.submit(debug=1, script_path='./test/core/test_file/test2.cmd')

def test_ClusterJob_preset():
    job = ClusterJob.config_job(
        commands = command_2_run,
        cluster = cluster,
        env_settings = [cluster.AMBER_GPU_ENV, cluster.G16_CPU_ENV],
        res_keywords = cluster.CPU_RES
    )
    print(job.sub_script_str)

@pytest.mark.accre # only run on ACCRE
def test_submit_job_id_ACCRE():
    '''
    only run on accre
    '''
    job = ClusterJob(accre.Accre(), sub_script_str=sub_script_str)
    job.submit(script_path='/home/shaoq1/EnzyHTP-test/test_job_manager/test.cmd')
    assert len(job.job_id) > 0
