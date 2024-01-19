"""Submodule contains code for charactization of the internal electric field of the enzyme.
+ field_strength_at()
    field strength at and along a specified bond. 

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


def field_strength_at(
        stru: Structure,
        p1: Union[ArrayLike, Atom],
        p2: Union[ArrayLike, Atom],
        d1: ArrayLike,
        location: str = 'center',
        region: str = None) -> float:
    """calculate the electric field strength of the internel electric field of an enzyme at
    a target bond or 2 points or 1 point and a direction (along the bond, at the {location}).
    Performs basic checks that the supplied Structure has charges. Output units are MV/cm.

    Args:
        stru:
            the structure of the target enzyme. MUST have charges initialized.
        p1, p2, d1:
            define the direction of the E projection
            and the point of measuring E.
            only 2 out of the 3 is needed. Allowed combinations are: (p1, p2), (p1, d1)
            p1, p2 can either be an Atom() or a point in the cartesian space.
        location:
            the location of the measurement when 2 points are specified.
            when it is (p1, d1), the p1 will be the position
            supported keywords: [center, p1] TODO add dipole_center
        region:
            the region in the structure that is considered as field source charges.
            use pymol selection syntax.

    Returns:
        The specified field strength in MV/cm."""

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

def ele_stab_energy() -> float:
    """"""
