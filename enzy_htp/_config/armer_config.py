"""Defines ARMerConfig() which holds configuration settings involved in ARMer process of interfacing
with high-performance clusters (HPCs). Controls default settings of different calculation types, etc.
Values can be set/accessed via string keys. A global instance of ARMerConfig() is owned by the module
singleton config variable. File also contains default_armer_config() which creates a default version of 
the ARMerConfig() object.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>

Date: 2023-09-22"""

from typing import Any, List
from copy import deepcopy

from .base_config import BaseConfig
from enzy_htp.core import _LOGGER

class ARMerConfig(BaseConfig):
    """Class that holds ARMer settings for enzy_htp. Similar to other classes in this
    sub-module, SHOULD NOT be directly created by the end users. Instead, it should be 
    accessed via the singleton config variable.
    Contain default settings of different job type that will be handled in ARMer. Modules
    such as /sampling or /energy should register the default values here.

    Attributes:
        QMCLUSTER_CPU_RES: dict() default resource settings for QM cluster jobs.
    """

    QMCLUSTER_CPU_RES: dict = {
        'core_type' : 'cpu',
        'nodes':'1',
        'node_cores' : '8',
        'job_name' : 'QMCluster_EnzyHTP',
        'partition' : '<fillthis>',
        'mem_per_core' : '3G', # in GB. Will used in gjf as downward round up after * 1024 and will - 1000 after multiple with cores
        'walltime' : '3-00:00:00',
        'account' : '<fillthis>',
        }
    """Default resource settings for QM cluster jobs from `energy`"""

    MD_CPU_RES: dict = {
        'core_type' : 'cpu',
        'nodes':'1',
        'node_cores' : '16',
        'job_name' : 'MD_EnzyHTP',
        'partition' : '<fillthis>',
        'mem_per_core' : '3G', # in GB
        'walltime' : '1-00:00:00',
        'account' : '<fillthis>'}
    """Default resource settings for CPU MD jobs from `geometry`"""

    MD_GPU_RES: dict = {
        'core_type' : 'gpu',
        'nodes':'1',
        'node_cores' : '1',
        'job_name' : 'MD_EnzyHTP',
        'partition' : '<fillthis>',
        'mem_per_core' : '8G', # in GB
        'walltime' : '3-00:00:00',
        'account' : '<fillthis>'}
    """Default resource settings for GPU MD jobs from `geometry`"""

    def required_executables(self) -> List[str]:
        """A list of all required executables for ARMerConfig."""
        return list()

    def required_env_vars(self) -> List[str]:
        """A list of all required environment variables for ARMerConfig."""
        return list()

    def required_py_modules(self) -> List[str]:
        """A list of all required environment python modules for ARMerConfig."""
        return list()


def default_armer_config() -> ARMerConfig:
    """Creates a deep-copied default version of the ARMerConfig() class."""
    return deepcopy(ARMerConfig())

