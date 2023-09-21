"""
"""
from typing import List, Any

from copy import deepcopy

from .base_config import BaseConfig


class AlphaFillConfig(BaseConfig):

    ALPHAFILL_EXE: str = "alphafill"

    CONFIG_FILE: str = "alphafill.conf"

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for AlphaFill."""
        return [self.ALPHAFILL_EXE]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for AlphaFill."""
        return []

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required python modules for AlphaFill."""
        return list()


def default_alphafill_config() -> AlphaFillConfig:
    """Creates a deep-copied default version of the AlphaFillConfig() class."""
    return deepcopy(AlphaFillConfig())
