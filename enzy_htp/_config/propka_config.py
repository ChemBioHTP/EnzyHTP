"""Submodule contains all configuration data for PROPKA.

Author: QZ Shao <shaoqz@icloud.com>

Date: 2025-08-07
"""
from copy import deepcopy
from typing import List

from .base_config import BaseConfig


class PropkaConfig(BaseConfig):
    """Holds all configuration data for PROPKA."""

    def required_env_vars(self) -> List[str]:
        return []

    def required_executables(self) -> List[str]:
        return []

    def required_py_modules(self) -> List[str]:
        return ['propka']

def default_propka_config() -> PropkaConfig:
    """Returns a deepcopy of the default PropkaConfig."""
    return deepcopy(PropkaConfig())

