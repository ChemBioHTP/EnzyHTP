"""Here is everything job manager need to know about ACCRE (Advanced Computing Center
for Research and Education) of Vanderbilt

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
import re
import os
from subprocess import CompletedProcess, SubprocessError, run
import time

from Class_Conf import Config
from ._interface import ClusterInterface


class Accre(ClusterInterface):
    '''
    The ACCRE interface
    '''
    #############################
    ### External use constant ###
    #############################
    NAME = 'ACCRE'

    # environment presets #
    AMBER_ENV = { 
        'CPU': '''module load GCC/6.4.0-2.28  OpenMPI/2.1.1
module load Amber/17-Python-2.7.14''', # only this version have sander.MPI
        'GPU': '''source /home/shaoq1/bin/amber_env/amber-accre.sh'''
    }

    G16_ENV = {
        'CPU':{ 'head' : '''module load Gaussian/16.B.01
mkdir $TMPDIR/$SLURM_JOB_ID
export GAUSS_SCRDIR=$TMPDIR/$SLURM_JOB_ID''',
                'tail' : '''rm -rf $TMPDIR/$SLURM_JOB_ID'''},
        'GPU': None
    }

    #############################
    ### Internal use constant ###
    #############################
    # command for submission
    SUBMIT_CMD = 'sbatch'
    # regex pattern for extracting job id from stdout
    JOB_ID_PATTERN = r'Submitted batch job ([0-9]+)'
    # command for control & monitor
    KILL_CMD = 'scancel'
    HOLD_CMD = 'scontrol hold'
    RELEASE_CMD = 'scontrol release'
    INFO_CMD = ['squeue', 'sacct'] # will check by order if previous one has no info
    # dict of job state
    JOB_STATE_MAP = {
        'pend' : ['CONFIGURING', 'PENDING', 'REQUEUE_FED', 'REQUEUE_HOLD', 'REQUEUED'],
        'run' : ['COMPLETING', 'RUNNING', 'STAGE_OUT'],
        'cancel' : ['CANCELLED', 'DEADLINE', 'TIMEOUT'],
        'complete' : ['COMPLETED'],
        'error' : ['BOOT_FAIL', 'FAILED', 'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 'REVOKED', 'STOPPED', 'SUSPENDED'],
        'exception': ['RESIZING', 'SIGNALING', 'SPECIAL_EXIT', 'RESV_DEL_HOLD' ]
    }

    RES_KEYWORDS_MAP = { 
        'core_type' : None,
        'nodes':'nodes=',
        'node_cores' : {'cpu': 'tasks-per-node=', 'gpu': 'gres=gpu:'},
        'job_name' : 'job-name=',
        'partition' : 'partition=',
        'mem_per_core' : {'cpu': 'mem-per-cpu=', 'gpu': 'mem-per-gpu='},
        'walltime' : 'time=',
        'account' : 'account='
    }

    PARTITION_VALUE_MAP = {
        'production' : {'cpu': 'production',
                        'gpu': '{maxwell,pascal,turing}'},
        'debug': {'cpu': 'debug',
                  'gpu': 'maxwell'}
    } # TODO do we really want to make partition general and parse it for each cluster?

    ##########################
    ### Submission Related ###
    ##########################
    @classmethod
    def parser_resource_str(cls, res_dict: dict) -> str:
        '''
        1. parser general resource keywords to accre specified keywords
        2. format the head of the submission script
        res_dict: the dictionary with general keywords and value
           (Available keys & value format:
                'core_type' : 'cpu',
                'node_cores' : '24',
                'job_name' : 'job_name',
                'partition' : 'production',
                'mem_per_core' : '4G',
                'walltime' : '24:00:00',
                'account' : 'xxx')
        return the string of the resource section
        '''
        parsered_res_dict = cls._parser_res_dict(res_dict)
        res_str = cls._format_res_str(parsered_res_dict)
        return res_str

    @classmethod
    def _parser_res_dict(cls, res_dict):
        '''
        parser keywords into ACCRE style
        '''
        new_dict = {}
        for k, v in res_dict.items():
            # remove core type line
            if k == 'core_type':
                continue
            # select core type
            if k in ('node_cores', 'mem_per_core'):
                new_k = cls.RES_KEYWORDS_MAP[k][res_dict['core_type']]
            else:
                new_k = cls.RES_KEYWORDS_MAP[k]
            new_dict[new_k] = v
        return new_dict

    @staticmethod
    def _format_res_str(parsered_res_dict):
        '''
        format parsered dictionary
        '''
        res_str = '#!/bin/bash\n'
        for k, v in parsered_res_dict.items():
            res_line = f'#SBATCH --{k}{v}\n'
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
        cmd = cls._format_submit_cmd(os.path.abspath(script_path))
        # debug
        if debug:
            print(cmd)
            return (cmd, sub_dir, script_path), None

        cwd = os.getcwd()
        # cd to sub_path
        os.chdir(sub_dir)
        try:    
            submit_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        except SubprocessError as e: # capture both error and timeout
            print(f'Error running {cmd}: {repr(e)}')
            print(f'stderr: {e.stderr}')
            print(f'stdout: {e.stdout}')
            os.chdir(cwd) # avoid messing up the dir
            raise e
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
    def get_job_info(cls, job_id: str, field: str, wait_time=3) -> str:
        '''
        get information about the job_id job by field keyword
        1. use squeue frist (fast) and
        if nothing is found (job finished)
        wait for wait time and 2. use sacct (slow)
        Arg:
            job_id
            field: supported keywords can be found at https://slurm.schedmd.com/sacct.html
            wait_time: for sacct run in second (default: 3s)
            **The `sacct` command takes some time (1-5s) to update the information of the job**   
        '''
        # get info
        # squeue
        cmd = f'{cls.INFO_CMD[0]} -u $USER -O JobID,{field}' # donot use the -j method to be more stable
        info_run = run(cmd, timeout=120, check=True,  text=True, shell=True, capture_output=True)
        # if exist
        info_out_lines = info_run.stdout.strip().splitlines()
        for info_line in info_out_lines:
            if job_id in info_line:
                job_field_info = info_line.strip().split()[1].strip().strip('+')
                return job_field_info
        # use sacct if squeue do not have info
        # wait a update gap
        if Config.debug > 1:
            print('No info from squeue. Switch to sacct')
        time.sleep(wait_time)
        cmd = f'{cls.INFO_CMD[1]} -j {job_id} -o {field}'
        try:
            info_run = run(cmd, timeout=60, check=True,  text=True, shell=True, capture_output=True)
        except SubprocessError as e:
            print(f'Error running {cmd}: {repr(e)}')
            print(f'stderr: {e.stderr}')
            print(f'stdout: {e.stdout}')
            raise e
            
        # if exist
        info_out = info_run.stdout.strip().splitlines()
        if len(info_out) >= 3:
            job_field_info = info_out[2].strip().strip('+')
            return job_field_info
        raise Exception(f'No information is found for {job_id}')
    
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
        raise Exception(f'Do not regonize state: {state}')