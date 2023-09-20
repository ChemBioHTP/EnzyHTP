"""This module defines the interface of MolDynParameterizer as an abstract class.
Concrete classes of MolDynParameterizer may show up in interfacing module of different
MD software.
TODO could it be a concrete class too that handle some general operations?

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-19
"""

from abc import ABC, abstractmethod

class MolDynParameterizer(ABC):
    """The parameterizer for Molecular Dynamics simulation."""
