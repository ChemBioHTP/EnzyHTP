"""Submodule contains code for charactization of the internal electric field of the enzyme.
+ field_strength_at()
    field strength at and along a specified bond.
+ ele_stab_energy()
    the electrostatic stablization energy. (See https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01424)

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: QZ Shao <shaoqz@icloud.com>

Date: 2022-11-06
"""
from typing import Union, List

from numpy.typing import ArrayLike
import numpy as np

from enzy_htp.core import _LOGGER
from enzy_htp.chemical import electric_field_strength
from enzy_htp.structure import Structure, Atom
from enzy_htp.structure.structure_operation.charge import init_charge
from enzy_htp.structure.structure_selection import select_stru


def ele_field_strength_at(
        stru: Structure,
        p1: Union[ArrayLike, Atom],
        p2: Union[ArrayLike, Atom] = None,
        d1: ArrayLike = None,
        location: str = 'center',
        region_pattern: str = None) -> float:
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
            supported keywords: [center, p1, p2] TODO support dipole_center
        region:
            the region in the structure that is considered as field source charges.
            use pymol selection syntax.

    Returns:
        The specified field strength in MV/cm."""
    # san check
    allowed_locs: List[str] = 'p1 p2 center'.split()
    if location not in allowed_locs:
        _LOGGER.error(f"The supplied location {location} is valid. Acceptable values include: '{', '.join(allowed_locs)}'")
        raise ValueError
    if not stru.has_charges():
        _LOGGER.info("The supplied structure DOES NOT have charges. Initializing...")
        init_charge(stru, ligand_fix="skip") #TODO finish ligand
    if p2 is None and d1 is None:
        raise Exception('Please provide a 2nd atom or point (p2=...) or a direction (d1=...)')

    # initialize variables
    stru_sele = select_stru(stru, region_pattern)
    if isinstance(p1, Atom):
        p1 = np.array(p1.coord, dtype=float)
    if p2 is not None:
        if isinstance(p2, Atom):
            p2 = np.array(p2.coord, dtype=float)
        d1 = p2 - p1
    d1 /= np.linalg.norm(d1)

    if location == 'center':
        p1 = (p1 + p2) / 2
    elif location == 'p1':
        p1 = p1
    elif location == 'p2':
        p1 = p2

    # calculate result
    result: float = 0.0
    for atom in stru_sele.atoms:
        source_coord = np.array(atom.coord)
        source_charge = atom.charge

        if source_charge is None:
            _LOGGER.error(f"{atom} dont have charge.")
            raise Exception

        result += electric_field_strength(source_coord, source_charge, p1, d1)

    return result

def ele_stab_energy() -> float:
    """"""
