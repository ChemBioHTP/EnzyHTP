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
    NAME = 'ACCRE'
    ##########################
    ### Submission Related ###
    ##########################
    # command for submission
    SUBMIT_CMD = 'sbatch'
    # regex pattern for extracting job id from stdout
    JOB_ID_PATTERN = r'Submitted batch job ([0-9]+)'
    # command for kill
    KILL_CMD = 'scancel'
    HOLD_CMD = 'scontrol hold'
    RELEASE_CMD = 'scontrol release'
    INFO_CMD = 'sacct'

    # dict of job state
    JOB_STATE_MAP = {
        'pend' : ['CONFIGURING', 'PENDING', 'REQUEUE_FED', 'REQUEUE_HOLD', 'REQUEUED'],
        'run' : ['COMPLETING', 'RUNNING', 'STAGE_OUT'],
        'cancel' : ['CANCELLED', 'DEADLINE', 'TIMEOUT'],
        'complete' : ['COMPLETED'],
        'error' : ['BOOT_FAIL', 'FAILED', 'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 'REVOKED', 'STOPPED', 'SUSPENDED'],
        'exception': ['RESIZING', 'SIGNALING', 'SPECIAL_EXIT', 'RESV_DEL_HOLD' ]
    }

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
    
    @classmethod
    def submit_job(cls, sub_dir, script_path, debug=0) -> tuple[str, str]:
        '''
        submit job submission script to the cluster queue
        cd to the sub_path and run submit cmd
        Return:
            (job_id, slurm_log_file_path)

        ACCRE sbatch rule:
            stdout: Submitted batch job ########
            file: slurm-#######.out will be generated in the *submission dir*
            exec: commands in the script will run under the *submission dir*
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
        slurm_log_path = cls._get_log_from_id(sub_dir, job_id)
        return (job_id, slurm_log_path)

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
        '''
        job_id = re.search(cls.JOB_ID_PATTERN, submit_job.stdout).group(1)
        return job_id

    @classmethod
    def _get_log_from_id(cls, sub_dir: str, job_id: str) -> str:
        '''
        file: slurm-#######.out will be generated in the *submission dir*
        '''
        file_path = sub_dir + f'/slurm-{job_id}.out'
        return file_path

    ###############################
    ### Post-submission Related ###
    ###############################

    @classmethod
    def kill_job(cls, job_id: str) -> CompletedProcess:
        cmd = f'{cls.KILL_CMD} {job_id}'
        kill_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        return kill_cmd

    @classmethod
    def hold_job(cls, job_id: str) -> CompletedProcess:
        cmd = f'{cls.HOLD_CMD} {job_id}'
        hold_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        return hold_cmd

    @classmethod
    def release_job(cls, job_id: str) -> CompletedProcess:
        cmd = f'{cls.RELEASE_CMD} {job_id}'
        release_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        return release_cmd

    @classmethod
    def get_job_info(cls, job_id: str, field: str) -> str:
        '''
        get information about the job_id job by field keyword
        supported keywords can be found at https://slurm.schedmd.com/sacct.html
        '''
        cmd = f'{cls.INFO_CMD} -j {job_id} -o {field}'
        info_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        job_field_info = info_cmd.stdout.strip().splitline()[2].strip()
        return job_field_info

    @classmethod
    def get_job_state(cls, job_id: str) -> tuple[str, str]:
        '''
        determine if the job is:
        Pend or Run or Complete or Canel or Error
        Return: 
            a tuple of
            (a str of pend or run or complete or canel or error,
                the real keyword form the cluster)
        '''
        state = cls.get_job_info(job_id, 'State')
        for k, v in cls.JOB_STATE_MAP.items():
            if state in v:
                return (k, state)