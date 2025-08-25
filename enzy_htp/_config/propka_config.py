"""Submodule contains all configuration data for PROPKA.

PROPKA is an empirical tool that predicts the pKa values of ionizable residues in proteins based on their 3D structure, aiding in the analysis of enzyme catalysis, protein stability, and charge states.

Author: QZ Shao <shaoqz@icloud.com>

Date: 2025-08-07
"""
from copy import deepcopy
from typing import List

from .base_config import BaseConfig


class PropkaConfig(BaseConfig):
    """Holds all configuration data for PROPKA.

    PROPKA computes the pKa values of ionizable amino acid residues in proteins using an empirical model, which is essential for accurate protonation state assignment in enzymatic and structural studies.
    """

    def required_env_vars(self) -> List[str]:
        return []

    def required_executables(self) -> List[str]:
        return []

    def required_py_modules(self) -> List[str]:
        return ['propka']

def default_propka_config() -> PropkaConfig:
    """Returns a deepcopy of the default PropkaConfig."""
    return deepcopy(PropkaConfig())

