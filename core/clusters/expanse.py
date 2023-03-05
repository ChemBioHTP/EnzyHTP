"""Here is everything job manager need to know about EXPANSE of SDSC
TODO: group them into a Slurm parent class

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
import re
import os
from subprocess import CompletedProcess, SubprocessError, run
import time

from Class_Conf import Config
from helper import round_by, run_cmd
from .accre import Accre


class Expanse(Accre):
    '''
    The EXPANSE interface
    '''
    #############################
    ### External use constant ###
    #############################
    NAME = 'EXPANSE'

    # environment presets #
    AMBER_ENV = { 
        'CPU': '''module load cpu/0.15.4  gcc/10.2.0  mvapich2/2.3.6
module load amber/20.21''', # only this version have sander.MPI
        'GPU': '''module load gpu/0.15.4 openmpi/4.0.4
module load amber/20'''
    }

    G16_ENV = {
        'CPU':{ 'head' : '''module load cpu/0.15.4
module load gaussian/16.C.01
export TMPDIR=/scratch/$USER/job_$SLURM_JOB_ID
mkdir $TMPDIR
export GAUSS_SCRDIR=$TMPDIR''',
                'tail' : '''rm -rf $TMPDIR'''},
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
        'node_cores' : {'cpu': 'ntasks-per-node=', 'gpu': 'gpus='},
        'job_name' : 'job-name=',
        'partition' : 'partition=',
        'mem_per_core' : {'cpu': 'mem=', 'gpu': 'mem='}, # previously using mem-per-gpu= change to mem= (calculate the total memory) base on issue #57
        'walltime' : 'time=',
        'account' : 'account='
    }

    PARTITION_VALUE_MAP = {
        'production' : {'cpu': 'shared',
                        'gpu': 'gpu-shared'},
        'debug': {'cpu': 'debug',
                  'gpu': 'gpu-debug'}
    } # TODO do we really want to make partition general and parse it for each cluster?

    @staticmethod
    def _format_res_str(parsered_res_dict):
        '''
        redefine the way to form res_str for expanse
        - add no-requeue
            SLURM will requeue jobs if there is a node failure.
            However, in some cases this might be detrimental if files get overwritten.
            So we use no-requeue here
        '''
        res_str = '#!/bin/bash\n'
        for k, v in parsered_res_dict.items():
            res_line = f'#SBATCH --{k}{v}\n'
            res_str += res_line
        res_str += '#SBATCH --no-requeue\n'
        return res_str

    @classmethod
    def _parser_res_dict(cls, res_dict):
        '''
        parser keywords into Expanse style
        '''
        new_dict = {}
        for k, v in res_dict.items():
            # remove core type line
            if k == 'core_type':
                if v == "gpu":
                    new_k = "ntasks-per-node="
                    new_v = 1 # default cpu for the gpu
                    new_dict[new_k] = new_v
                continue
            # select core type
            if k in ('node_cores', 'mem_per_core'):
                new_k = cls.RES_KEYWORDS_MAP[k][res_dict['core_type']]
            else:
                new_k = cls.RES_KEYWORDS_MAP[k]
            new_dict[new_k] = v
        # process total mem
        for k in new_dict:
            if k == 'mem=':
                mem_per_core_n_gb = new_dict[k].rstrip('GB')
                total_mem = round_by(float(mem_per_core_n_gb) * float(res_dict['node_cores']), 0.1) # round up
                new_dict[k] = f'{total_mem}G'
        return new_dict

