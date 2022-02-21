# python imports
import sys
import logging

# enzy_htp imports
from .structure import Structure, Residue
from .core import EnvironmentManager, MissingEnvironmentElement, InvalidResidueCode
from .chemical import ResidueType


def init_logging() -> None:
    logging.basicConfig(
        level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)]
    )


em = EnvironmentManager(
    env_vars=["$AMBERHOME"],
    executables=[
        "tleap",
        "ambpdb",
        "cpptraj",
        "mpirun",
        "$AMBERHOME/bin/sander.MPI",
        "$AMBERHOME/bin/pmemd.cuda",
        "$AMBERHOME/bin/MMPBSA.py.MPI",
    ],
)

init_logging()
em.check_environment()
