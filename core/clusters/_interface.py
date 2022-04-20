"""Interface regulation of supporting different clusters

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
from abc import ABC, abstractmethod

class ClusterInterface(ABC):
    '''
    Defines the interface of a cluster:
    property
    '''
    @property
    @abstractmethod
    def SUBMIT_CMD(self) -> str:
        pass