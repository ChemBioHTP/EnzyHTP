"""Manage the job queue for running jobs on a cluster. Interface with different linux resource managers (e.g. slurm)
The general workflow is
    (take job commands)
    1) compile the submission script 
    2) submit and monitor the job 
    3) record the completion of the job.
    (give files containing the stdout/stderr)
In a dataflow point of view it should be analogous to subprocess.run()

Feature:
    - Allow users to add support for their own clusters. (By making new ClusterInterface classes)

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
from typing import Union
from subprocess import run

from .clusters import *  # pylint: disable=unused-wildcard-import,wildcard-import
from core.clusters._interface import ClusterInterface

class ClusterJob():
    '''
    This class handle jobs for cluster calculation
    API:
    property:
        # command related #
        commands:   a list of commands to run
        env_modules:modules that set environment on the cluster (set by module)
        env_lines:  commnad lines that set environment on the cluster (set by line)
        # cluster related #
        cluster:    cluster used for running the job (pick from list in /core/cluster/)
        res_keyword: (set by keyword)
            core_type:  gpu | cpu (TODO(shaoqz): support gpu&cpu)
            partation:  a list supporting `core_type` specified by each cluster
            account:    account for using `partation`
            nodes:      number nodes
            node_cores: number of cores of each node
            mem_per_core: memory per core
            walltime:   time limit of the job
        res_lines:  a cluster specified line (set resources by line)
        sub_script: directly the submission script
    method:

    '''

    def __init__(self,
                 commands: list[str] = None,
                 env_modules: list[str] = None,
                 env_lines: list[str] = None,
                 cluster: Union[ClusterInterface, None] = None,
                 res_keyword: dict[str, str] = None,
                 res_lines: list[str] = None,
                 sub_script: str = None
                 ) -> None:
        self.commands = commands
        self.env_modules = env_modules
        self.env_lines = env_lines
        self.cluster = cluster
        self.res_keyword = res_keyword
        self.res_lines = res_lines
        self.sub_script = sub_script
        
    ### method ###
    def submit(self, debug=0):
        '''
        submit current job to the cluster queue
        '''
        cmd = self.cluster.submit_cmd

        if not debug:
            run(cmd, check=True, text=True, shell=True, capture_output=True)
        else:
            print(cmd)

    def config_job(self):
        '''
        config job with:
        1. resource settings
        2. environment
        3. commands to run
        '''
        pass

    def make_sub_script(self):
        pass

