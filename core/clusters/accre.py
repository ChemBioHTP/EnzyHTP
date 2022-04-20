"""Here is everything job manager need to know about ACCRE (Advanced Computing Center
for Research and Education) of Vanderbilt

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
from ._interface import ClusterInterface


class Accre(ClusterInterface):
    '''
    The ACCRE interface
    API:
    SUBMIT_CMD      # job script submission command

    '''

    SUBMIT_CMD = 'sbatch'

    # common environment setting
    AMBER_GPU_ENV = [
        'export AMBERHOME=/dors/csb/apps/amber19/', 'export CUDA_HOME=$AMBERHOME/cuda/10.0.130',
        'export LD_LIBRARY_PATH=$AMBERHOME/cuda/10.0.130/lib64:$AMBERHOME/cuda/RHEL7/10.0.130/lib:$LD_LIBRARY_PATH'
    ]

    # route format
    # support for gaussian
    # support for amber
    def __init__(self) -> None:
        pass
