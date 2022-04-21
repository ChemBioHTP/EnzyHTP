"""Interface regulation of supporting different clusters

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
from abc import ABC, abstractmethod
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
    def SUBMIT_CMD(self) -> str:
        '''
        cmd for submission TODO dont need this
        '''
        pass

    ### classmethods ###
    @classmethod
    @abstractmethod
    def format_resource_str(cls, res_dict: dict) -> str:
        '''
        format the head of the submission script
        res_dict: the dictionary with exact the keyword and value
                  required by the cluster
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
