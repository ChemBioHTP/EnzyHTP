"""Submodule contains code for calculating the electric field strength in an enzyme system along a specified
axis. Calculates field strength for an enzyme using Amber RESP charges. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: QZ Shao <shaoqz@icloud.com>

Date: 2022-11-06
"""
from typing import Union, List

from numpy.typing import ArrayLike
import numpy as np
import pandas as pd

from enzy_htp.core import _LOGGER
from enzy_htp.structure import Structure, Atom


def electric_field(stru: Structure,
                   p1: Union[ArrayLike, Atom],
                   p2: Union[ArrayLike, Atom],
                   location: str = 'center',
                   mask: str = None) -> float:
    """Top level method for calculating the electric field strength of an enzyme, given an atom mask, two points or atoms,
    and the location between the two points or atoms. Note that only the stru and p1/p2 are required. Performs basic checks
    that the supplied Structure has charges. Output units are MV/cm.

    Args:
        stru: The Structure to use. MUST have charges.
        p1: The point or Atom to use as the first point for the direction vector.
        p2: The point or atom to use as the second point for the direction vector.
        location: Location to take the electric field strength at along the direction vector. Allowed values are p1, p2, and center.
        mask: The selection mask to use. If not supplied, all atoms are used.

    Returns:
        The specified field strength in MV/cm.

    """

    allowed_locs: List[str] = 'p1 p2 center'.split()
    if location not in allowed_locs:
        _LOGGER.error(f"The supplied location {location} is valid. Acceptable values include: '{', '.join(allowed_locs)}'")
        exit(1)

    if not stru.has_charges():
        _LOGGER.error("The supplied structure DOES NOT have charges. Exiting...")
        exit(1)
        pass

    if type(p1) == Atom:
        p1 = p1.coord()

    if type(p2) == Atom:
        p2 = p2.coord()

    p1 = p1.astype(float)
    p2 = p2.astype(float)

    direction = p2 - p1

    direction /= np.linalg.norm(direction)

    if location == 'center':
        p1 = (p1 + p2) / 2
    elif location == 'p1':
        p1 = p1
    elif location == 'p2':
        p1 = p2

    result: float = 0.0
    k: float = 332.4

    for atom in stru.atoms:
        #TODO(CJ): add in the masking of the structure
        p0: np.ndarray = np.array(atom.coord)
        r: np.ndarray = p1 - p0
        r_m: float = np.linalg.norm(r)
        r_u: np.ndarray = r / r_m
        result += np.dot(((k * atom.charge) / (r_m**2)) * r_u, direction)

    return result
