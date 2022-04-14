"""Here is everything job manager need to know about ACCRE (Advanced Computing Center
for Research and Education) of Vanderbilt

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
from ._interface import ClusterInterface

class Accre(ClusterInterface):
    '''
    The ACCRE interface
    '''
    # route format
    # support for gaussian
    # support for amber
    def __init__(self) -> None:
        pass

    submit_cmd = 'sbatch'
    # def submit_cmd(self):
    #     '''
    #     return the cmd used for submission
    #     '''
    #     return ''