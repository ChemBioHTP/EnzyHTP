"""TODO(CJ)"""


from typing import List, Any

from copy import deepcopy

from .base_config import BaseConfig


class Mole2Config(BaseConfig):



    PROBE:float=1.0

    INNER:float=1.0

    MOLE2:float=""


    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Mole2."""
        return [self.MOLE2]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for Mole2."""
        return []

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required python modules for Mole2."""
        return list()


def default_mole2_config() -> Mole2Config:
    """Creates a deep-copied default version of the Mole2Config() class."""
    return deepcopy(Mole2Config())
