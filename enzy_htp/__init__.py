# python imports
import sys
import logging
from typing import Tuple

# enzy_htp imports
from .config import Config

# from .interface import Interface
from .structure import Structure, Residue
from .core import (
    EnvironmentManager,
    MissingEnvironmentElement,
    InvalidResidueCode,
    _LOGGER,
)
from .preparation import PDBLine, PDBPrepper, read_pdb_lines
from .chemical import ResidueType

from .molecular_mechanics import AmberInterface

CONFIG = None
INTERFACE = None


# def init() -> Tuple[Config, Interface]:
#    """"""
#    global CONFIG
#    global INTERFACE
#    CONFIG = Config()
#    #INTERFACE = Interface(CONFIG)
#    return (CONFIG, INTERFACE)
