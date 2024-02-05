from typing import List, Any

from copy import deepcopy

from .base_config import BaseConfig


class ModellerConfig(BaseConfig):



    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for AlphaFill."""
        return []

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for AlphaFill."""
        return []

    def required_py_modules(self) -> List[str]:
        """A hardcoded list of required python modules for AlphaFill."""
        return ["modeller"]

def default_modeller_config() -> ModellerConfig:
    """Creates a deep-copied default version of the ModellerConfig() class."""

    return deepcopy(ModellerConfig())
