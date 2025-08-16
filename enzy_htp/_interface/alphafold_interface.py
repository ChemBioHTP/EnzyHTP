"""Defines an AlphafoldInterface class that serves as a bridge for enzy_htp to utilize Alphafold software.

Author: Gemini
Date: 2025-08-16
"""
from __future__ import annotations

from .base_interface import BaseInterface
from enzy_htp.structure import Structure
from enzy_htp._config.alphafold_config import AlphafoldConfig

class AlphafoldInterface(BaseInterface):
    """Class that provides a direct inteface for enzy_htp to utilize Alphafold software.
    """

    def __init__(self, parent, config: AlphafoldConfig = None) -> None:
        """Simplistic constructor that optionally takes an AlphafoldConfig object as its only argument.
        Calls parent class."""
        super().__init__(parent, config, AlphafoldConfig)

    def run(self, sequence: str, **kwargs) -> Structure:
        """Run Alphafold to predict a structure from a sequence."""
        raise NotImplementedError
