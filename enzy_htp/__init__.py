# python imports
import sys
import logging

# enzy_htp imports
from .config import CONFIG
from .interface import Interface
from .structure import Structure, Residue
from .core import (
    EnvironmentManager,
    MissingEnvironmentElement,
    InvalidResidueCode,
    _LOGGER,
)
from .preparation import PDBLine, PDBPrepper, read_pdb_lines
from .chemical import ResidueType

# def init_logging() -> None:
#    logging.basicConfig(
#        level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)], foj
#    )
# em.check_environment()

INTERFACE = Interface(CONFIG)
