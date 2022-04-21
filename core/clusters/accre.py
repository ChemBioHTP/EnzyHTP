"""Here is everything job manager need to know about ACCRE (Advanced Computing Center
for Research and Education) of Vanderbilt

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
from ._interface import ClusterInterface


class Accre(ClusterInterface):
    '''
    The ACCRE interface
    API:
    SUBMIT_CMD      # job script submission command

    '''
    # submission format
    SUBMIT_CMD = 'sbatch'

    # environment presets
    AMBER_GPU_ENV = '''export AMBERHOME=/dors/csb/apps/amber19/
export CUDA_HOME=$AMBERHOME/cuda/10.0.130
export LD_LIBRARY_PATH=$AMBERHOME/cuda/10.0.130/lib64:$AMBERHOME/cuda/RHEL7/10.0.130/lib:$LD_LIBRARY_PATH'''

    G16_CPU_ENV = '''module load Gaussian/16.B.01
mkdir $TMPDIR/$SLURM_JOB_ID
export GAUSS_SCRDIR=$TMPDIR/$SLURM_JOB_ID''' # remember to add the command that remove the SCRDIR

    # resource preset TODO
    CPU_RES = { 'core_type' : 'cpu',
                'node_cores' : '24',
                'job_name' : 'job_name',
                'partition' : 'production',
                'mem_per_core' : '4G',
                'walltime' : '24:00:00',
                'account' : 'xxx'}

    # resource format
    def format_resource_str(self, res_dict) -> str:
        '''
        format the head of the submission script for ACCRE
        res_dict: the dictionary with exact the keyword and value
                  required by ACCRE
        '''
        res_str = '#!/bin/bash\n'
        for k, v in res_dict.items():
            res_line = f'#SBATCH --{k}={v}\n'
            res_str += res_line
        return res_str