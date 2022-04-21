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
from typing import Union
from plum import dispatch

from .clusters import *  # pylint: disable=unused-wildcard-import,wildcard-import
from core.clusters._interface import ClusterInterface
from helper import line_feed
from Class_Conf import Config

class ClusterJob():
    '''
    This class handle jobs for cluster calculation
    API:
    constructor:
        ClusterJob.config_job()
    property:
        cluster:    cluster used for running the job (pick from list in /core/cluster/)
        sub_script_str: submission script content
        sub_script_path: submission script path
    method:
        submit()
    '''

    def __init__(self,
                 cluster: ClusterInterface,
                 sub_script_str: str
                 ) -> None:
        self.cluster = cluster
        self.sub_script_str = sub_script_str
        
        self.sub_script_path: str

    ### method ### 
    def submit(self, debug=0, script_path=None):
        '''
        submit the job to the cluster queue. Make the submission script. Submit.
        '''
        self._deploy_sub_script(script_path)
        self._submit(debug=debug)

    def _deploy_sub_script(self, out_path: str) -> None:
        '''
        deploy the submission scirpt for current job
        store the out_path to self.sub_script_path
        '''
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write(self.sub_script_str)
        self.sub_script_path = out_path
    
    def _submit(self, debug=0) -> None:
        '''
        interal method that submit current job submission script to the cluster queue
        '''
        cmd = ' '.join((self.cluster.SUBMIT_CMD, self.sub_script_path))

        if not debug:
            run(cmd, check=True, text=True, shell=True, capture_output=True)
        else:
            print(cmd)

    @classmethod
    def config_job( cls, 
                commands: Union[list[str], str],
                cluster: ClusterInterface,
                env_settings: Union[list[str], str],
                res_keywords: dict[str, str]
                ) -> 'ClusterJob':
        '''
        config job and generate a ClusterJob instance

        Args:
        commands: 
            commands to run. Can be a str of commands or a list containing strings of commands.
        cluster: 
            cluster for running the job. Should be a ClusterInterface object. 
            Available clusters can be found under core/clusters as python class defination.
            To define a new cluster class for support, reference the ClusterInterface requirement.
        env_settings: 
            environment settings in the submission script. Can be a string or list of strings
            for cmds in each line.
            Since environment settings are attached to job types. It is more conserved than the command.
            **Use presets in ClusterInterface classes to save effort**
        res_keywords: 
            resource settings. Can be a dictionary indicating each keywords or the string of the whole section.
            The name and value should be exactly the same as required by the cluster.
            **Use presets in ClusterInterface classes to save effort**
        
        Return:
        A ClusterJob object
        '''
        command_str = cls._get_command_str(commands)
        env_str = cls._get_env_str(env_settings)
        res_str = cls._get_res_str(res_keywords, cluster)
        sub_script_str = cls._get_sub_script_str(
                            command_str, 
                            env_str, 
                            res_str, 
                            f'# {Config.WATERMARK}{line_feed}'
                            )

        return cls(cluster, sub_script_str)

    @staticmethod
    @dispatch
    def _get_command_str(cmd: list) -> str:
        return line_feed.join(cmd)+line_feed

    @staticmethod
    @dispatch
    def _get_command_str(cmd: str) -> str:
        return cmd+line_feed

    @staticmethod
    @dispatch
    def _get_env_str(env: list) -> str:
        return line_feed.join(env)+line_feed

    @staticmethod
    @dispatch
    def _get_env_str(env: str) -> str:
        return env+line_feed

    @staticmethod
    @dispatch
    def _get_res_str(res: dict, 
                     cluster: ClusterInterface) -> str:
        return cluster.format_resource_str(res)

    @staticmethod
    @dispatch
    def _get_res_str(res: str, 
                     cluster: ClusterInterface) -> str:
        return res

    @staticmethod
    def _get_sub_script_str(command_str: str, env_str: str, res_str: str, watermark: str) -> str:
        '''
        combine command_str, env_str, res_str to sub_script_str
        '''
        sub_script_str= line_feed.join((res_str, watermark, env_str, command_str))
        return sub_script_str

    @dispatch
    def _(self):
        '''
        dummy method for dispatch
        '''
        pass