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
from subprocess import run
from plum import dispatch

from .clusters import *  # pylint: disable=unused-wildcard-import,wildcard-import
from core.clusters._interface import ClusterInterface
from helper import mkdir, line_feed

class ClusterJob():
    '''
    This class handle jobs for cluster calculation
    API:
    property:
        # # command related #
        # commands:   a list of commands to run
        # env_modules:modules that set environment on the cluster (set by module)
        # env_lines:  commnad lines that set environment on the cluster (set by line)
        # # cluster related #
        # cluster:    cluster used for running the job (pick from list in /core/cluster/)
        # res_keyword: (set by keyword)
        #     core_type:  gpu | cpu (TODO(shaoqz): support gpu&cpu)
        #     partation:  a list supporting `core_type` specified by each cluster
        #     account:    account for using `partation`
        #     nodes:      number nodes
        #     node_cores: number of cores of each node
        #     mem_per_core: memory per core
        #     walltime:   time limit of the job
        # res_lines:  a cluster specified line (set resources by line)
        # sub_script: directly the submission script
    method:
    ClusterJob():   Usage:
                    1. use sub_script_str to create a clusterjob. overwrite any other related configs OR
                    2. create a job first so that can config with useful cluster presets
    '''

    def __init__(self,
                 name: str,
                 cluster: ClusterInterface,
                 sub_script_str: str = None
                 ) -> None:
        self.name = name
        self.cluster = cluster
        self.sub_script_str = sub_script_str
        
        self.commands: str
        self.sub_script_path: str
        self.command_str: str
        self.env_str: str
        self.res_str: str
        
    ### method ### 
    def _deploy_sub_script(self, out_path: str) -> None:
        '''
        deploy the submission scirpt for current job
        store the out_path to self.sub_script_path
        '''
        self.get_sub_scrpit_str()
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write(self.sub_script_str)
        self.sub_script_path = out_path
    
    def get_sub_scrpit_str(self) -> str:
        '''
        get str for submission script.
        Generate with self.command_str, self.env_str, self.res_str if does not exist

        return the reference of the self.sub_script_str
        '''
        if self.sub_script_str is None:
            self.sub_script_str = line_feed.join((self.res_str, self.env_str, self.command_str))

        return self.sub_script_str

    def _submit(self, debug=0) -> None:
        '''
        interal method that submit current job submission script to the cluster queue
        '''
        cmd = ' '.join((self.cluster.SUBMIT_CMD, self.sub_script_path))

        if not debug:
            run(cmd, check=True, text=True, shell=True, capture_output=True)
        else:
            print(cmd)

    def submit(self, debug=0, script_path=None):
        '''
        submit the job to the cluster queue. Make the submission script. Submit.
        '''
        if script_path is None:
            mkdir('./tmp')
            script_path = './tmp/'+self.name+'.cmd'

        self._deploy_sub_script(script_path)
        self._submit(debug=debug)

    @dispatch
    def config( self, 
                command: list[str],
                env_setting: list[str],
                res_keyword: dict[str, str]
                ):
        '''
        config job with:
        1. resource settings
        2. environment
        3. commands to run
        '''
        self._get_command_str(command)

    @dispatch
    def _get_command_str(self, cmd_input: list[str]):
        pass

    @dispatch
    def _get_env_str(self, env_input: list[str]):
        pass

    @dispatch
    def _get_res_str(self, res_input: dict[str, str]):
        pass

