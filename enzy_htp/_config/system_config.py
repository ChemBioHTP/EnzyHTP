"""Defines SystemConfig() which holds configuration settings for processes that will be spawned by 
enzy_htp. Controls number of cores, memory usage, and locations of temp directories. Values can
be set/accessed via string keys. A global instance of SystemConfig() is owned by the module singleton
config variable. File also contains default_system_config() which creates a default version of 
the SystemConfig() object.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-10-16
"""
import os
from typing import Any
from copy import deepcopy

from enzy_htp.core import _LOGGER


class SystemConfig:
    """Class that holds system settings for enzy_htp. Similar to other classes in this
    sub-module, SHOULD NOT be directly created by the end users. Instead, it should be 
    accessed via the singleton config variable. 

    Attributes:
        N_CORES: int() describing max number of cores enzy_htp has access to.
        MEM_PER_CORE: int() max memory per core in megabytes (1GB = 1000).
        WORK_DIR: str() saying which directory enzy_htp is being run from.
        SCRATCH_DIR: str() with default temp directory. Set to WORK_DIR/scratch by default.
    """

    N_CORES: int = 24
    """Number of cores each process managed by enzy_htp has access to."""

    MEM_PER_CORE: int = 2000
    """Memory in megabytes that each enzy_htp process has."""

    WORK_DIR: str = os.getcwd()
    """Directory work is being done in by enzy_htp."""

    SCRATCH_DIR: str = f"{WORK_DIR}/scratch"
    """Default temp directory for enzy_htp. Defaults to WORK_DIR/scratch"""

    def __getitem__(self, key: str) -> Any:
        """Accessor for SystemConfig() that leverages [] operator syntax.
        
        Args:
            key: a str() key with the name of the variable you would like to access.

        Returns:
        """
        if key.count("."):
            _LOGGER.error(
                f"No nested variables currently exist in the SystemConfig() class. Exiting..."
            )
            exit(1)
        else:
            return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        """Setter for SystemConfig() that leverages [] operator syntax.

        Args:
            key: a str() with the name of the variable you would like to update.
            value: Whatever value you want to set the appropriate attribute to.

        Returns:
            Nothing
        """
        if key.count("."):
            _LOGGER.error(
                f"No nested variables currently exist in the SystemConfig() class. Exiting..."
            )
            exit(1)

        else:
            setattr(self, key, value)


def default_system_config() -> SystemConfig:
    """Creates a deep-copied default version of the SystemConfig() class."""
    return deepcopy(SystemConfig())
