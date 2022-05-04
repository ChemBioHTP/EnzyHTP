"""Interface regulation of supporting different clusters

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
from abc import ABC, abstractmethod
from ctypes import Union
from subprocess import CompletedProcess

class ClusterInterface(ABC):
    '''
    Defines the interface of a cluster
    ----------
    SUBMIT_CMD: the command for job submission

    '''
    ### class attribute ###
    @property
    @abstractmethod
    def NAME(self) -> str:
        pass

    @property
    @abstractmethod    
    def G16_ENV(self) -> dict:
        '''
        the cluster need to know the environment setting for each software
        This should be a dictionary with CPU and GPU as key.
        In each value, you can also use head and tail for env setting before and after the command.
        (see ClusterJob.config_job for more detail)
        '''
        pass

    @property
    @abstractmethod    
    def AMBER_ENV(self) -> dict:
        '''
        the cluster need to know the environment setting for each software
        This should be a dictionary with CPU and GPU as key.
        In each value, you can also use head and tail for env setting before and after the command.
        (see ClusterJob.config_job for more detail)
        '''
        pass
    
   ### classmethods ###
    @classmethod
    @abstractmethod
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
        pass

    @classmethod
    @abstractmethod
    def submit_job(cls, sub_dir: str, script_path: str) -> tuple[str, str]:
        '''
        submit job submission script to the cluster queue
        cd to the sub_path and run submit cmd
        return (job_id, cluster_log_file_path)
        '''
        pass

    @classmethod
    @abstractmethod
    def kill_job(cls, job_id: str) -> CompletedProcess:
        pass

    @classmethod
    @abstractmethod
    def hold_job(cls, job_id: str) -> CompletedProcess:
        pass

    @classmethod
    @abstractmethod
    def release_job(cls, job_id: str) -> CompletedProcess:
        pass

    @classmethod
    @abstractmethod
    def get_job_state(cls, job_id: str) -> tuple[str, str]:
        '''
        determine if the job is:
        Pend or Run or Complete or Canel or Error
        Return:
            a tuple of
            (a str of pend or run or complete or canel or error,
             the real keyword form the cluster)
        '''
        pass
   
