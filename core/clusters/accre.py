"""Here is everything job manager need to know about ACCRE (Advanced Computing Center
for Research and Education) of Vanderbilt

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
import re
import os
from subprocess import CompletedProcess, run
from ._interface import ClusterInterface


class Accre(ClusterInterface):
    '''
    The ACCRE interface
    '''
    ##########################
    ### Submission Related ###
    ##########################
    # command for submission
    SUBMIT_CMD = 'sbatch'

    ### presets ###
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

    @classmethod
    def format_resource_str(cls, res_dict: dict) -> str:
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
    
    ###############################
    ### Post-submission Related ###
    ##########################
    # regex pattern for extracting job id from stdout
    JOB_ID_PATTERN = r'Submitted batch job ([0-9]+)'

    @classmethod
    def submit_job(cls, sub_dir, script_path, debug=0) -> None:
        '''
        interal method that submit current job submission script to the cluster queue
        cd to the sub_path and run submit cmd
        '''
        cmd = cls._format_submit_cmd(script_path)
        # debug
        if debug:
            print(cmd)
            return cmd

        cwd = os.getcwd()
        # cd to sub_path
        os.chdir(sub_dir)        
        submit_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        # TODO(shaoqz) timeout condition is hard to test
        # cd back
        os.chdir(cwd)
        
        job_id = cls._get_job_id_from_submit(submit_cmd)
        return job_id

    @classmethod
    def _format_submit_cmd(cls, sub_script_path: str) -> str:
        '''
        format the command for job submission on ACCRE
        return the format command.
        '''
        cmd = ' '.join((cls.SUBMIT_CMD, sub_script_path))
        return cmd

    @classmethod
    def _get_job_id_from_submit(cls, submit_job: CompletedProcess) -> str:
        '''
        extract job id from the output of the submission command run.
        ACCRE sbatch rule:
            stdout: Submitted batch job ########
            file: slurm-#######.out will be generated in the *submission dir*
            exec: commands in the script will run under the *submission dir*
        '''
        job_id = re.search(cls.JOB_ID_PATTERN, submit_job.stdout).group(1)
        return job_id

    @classmethod
    def kill_job(cls, job_id: str) -> str:
        pass
