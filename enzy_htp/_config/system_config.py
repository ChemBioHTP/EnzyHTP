"""
TODO(CJ)
"""
from typing import Any
from copy import deepcopy

# TODO(CJ): need documentation here
class SystemConfig:

    N_CORES = 24

    MEM_PER_CORE = 2000

    pass

    def __getitem__(self, key: str) -> Any:
        if key.count("."):
            key1, key2 = key.split(".")
            return getattr(self, key1)[key2]
        else:
            return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        if key.count("."):
            key1, key2 = key.split(".")
            GaussianConfig.__dict__[key1][key2] = value
        else:
            setattr(self, key, value)


def default_system_config() -> SystemConfig:
    return deepcopy(SystemConfig())
