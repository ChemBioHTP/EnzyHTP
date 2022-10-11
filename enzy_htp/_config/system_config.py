"""
TODO(CJ)
"""
import os
from typing import Any
from copy import deepcopy

# TODO(CJ): need documentation here
class SystemConfig:

    N_CORES = 24

    MEM_PER_CORE = 2000

    WORK_DIR = os.getcwd()

    SCRATCH_DIR = f"{WORK_DIR}/scratch"

    def __getitem__(self, key: str) -> Any:
        if key.count("."):
            key1, key2 = key.split(".")
            return getattr(self, key1)[key2]
        else:
            return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        if key.count("."):
            key1, key2 = key.split(".")
            SystemConfig.__dict__[key1][key2] = value
        else:
            setattr(self, key, value)


def default_system_config() -> SystemConfig:
    return deepcopy(SystemConfig())
