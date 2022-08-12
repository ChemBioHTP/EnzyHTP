"""


TODO(CJ)
"""

from pathlib import Path
from typing import List, Dict

from enzy_htp import interface
from enzy_htp.structure import Structure
from enzy_htp.core import UnsupportedMethod


def _sample_amber_md(pdb: str, point: int) -> List[Structure]:
    """ """
    pass
    work_dir = Path(pdb).parent  # TODO(CJ): change this
    (prmtop, inpcrd) = interface.amber.build_param_files(pdb, work_dir)
    print(prmtop)
    print(inpcrd)
    fname = interface.amber.md_run(prmtop, inpcrd, work_dir)
    print(fname)


def sample_geometries(pdb: str, point: int, engine: str) -> List[Structure]:
    """

    Args:
        pdb: a string containing the path to a structure in .pdb file format.
        point:
        engine:

    Returns:
        A list() containing the sampled geometries.
    """
    IMPLEMENTATION: Dict = {"Amber_MD": _sample_amber_md}

    if engine in IMPLEMENTATION:
        return IMPLEMENTATION[engine](pdb, point)
    else:
        raise UnsupportedMethod(
            f"{engine} is not supported for sample_geometries(). Supported methods include: {', '.join(list(IMPLEMENTATION.keys()))}"
        )
