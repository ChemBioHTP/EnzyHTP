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
from scipy.constants import physical_constants

from numpy.typing import ArrayLike
import numpy as np

from enzy_htp.core import _LOGGER
from enzy_htp.chemical import electric_field_strength
from enzy_htp.structure import Structure, Atom
from enzy_htp.structure.structure_operation.charge import init_charge
from enzy_htp.structure.structure_selection import select_stru


def ele_field_strength_at_along(
        stru: Structure,
        p1: Union[ArrayLike, Atom], #TODO support atom key
        p2: Union[ArrayLike, Atom] = None,
        d1: ArrayLike = None,
        location: str = 'center',
        region_pattern: str = None,
        unit: str = "kcal/(mol*e*Ang)",) -> float:
    """At {location}, along {d1} or {p2 - p1}
    calculate the electric field strength of the internel electric field of an enzyme at
    a target bond or 2 points or 1 point and a direction (along the bond, at the {location}).
    Performs basic checks that the supplied Structure has charges.

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
        unit:
            the unit of the output EF

    Returns:
        The specified field strength in {unit}."""
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
        p1 = stru.get_corresponding_atom(p1)
        p1 = np.array(p1.coord, dtype=float)
    if p2 is not None:
        if isinstance(p2, Atom):
            p2 = stru.get_corresponding_atom(p2)
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

        result += electric_field_strength(source_coord, source_charge, p1, d1, unit=unit)

    return result

def ele_field_strength_at(
        stru: Structure,
        p1: Union[ArrayLike, Atom],
        p2: Union[ArrayLike, Atom] = None,
        location: str = None,
        region_pattern: str = None,
        unit: str = "kcal/(mol*e*Ang)",) -> ArrayLike:
    """At {location}
    calculate the electric field strength of the internel electric field of an enzyme at
    a target bond or 2 points or 1 point.
    Performs basic checks that the supplied Structure has charges.
    NOTE: this is a seperate function because we want to avoid overload that returns
    different data types. the caller should explicitly know what reture type to expect.

    Args:
        stru:
            the structure of the target enzyme. MUST have charges initialized.
        p1, p2, location:
            define the point of measuring E.
            By default use p1, if location is specifcied, p2 is needed.
            p1, p2 can either be an Atom() or a point in the cartesian space.
            supported location: [center] TODO support dipole_center
        region:
            the region in the structure that is considered as field source charges.
            use pymol selection syntax.
        unit:
            the unit of the output EF

    Returns:
        The specified field strength in {unit}."""
    # san check
    allowed_locs: List[str] = ["center"]
    if location not in allowed_locs:
        _LOGGER.error(f"The supplied location {location} is valid. Acceptable values include: '{', '.join(allowed_locs)}'")
        raise ValueError
    if not stru.has_charges():
        _LOGGER.info("The supplied structure DOES NOT have charges. Initializing...")
        init_charge(stru, ligand_fix="skip") #TODO finish ligand
    if location is not None and p1 is None:
        raise Exception('Please provide a 2nd atom when location is used')

    # initialize variables
    stru_sele = select_stru(stru, region_pattern)
    if isinstance(p1, Atom):
        p1 = np.array(p1.coord, dtype=float)
    if location == 'center':
        if isinstance(p2, Atom):
            p2 = np.array(p2.coord, dtype=float)
        p1 = (p1 + p2) / 2

    # calculate result
    result: float = np.array((0.0, 0.0, 0.0))
    for atom in stru_sele.atoms:
        source_coord = np.array(atom.coord)
        source_charge = atom.charge

        if source_charge is None:
            _LOGGER.error(f"{atom} dont have charge.")
            raise Exception

        result += electric_field_strength(source_coord, source_charge, p1, unit=unit)

    return result

def ele_stab_energy_of_bond(
        bond_dipole: float, 
        bond_field_strength: float,) -> float:
    """calculate the electrostatic stablization energy of
    a bond dipole. Assume the dipole is in the same or opposite
    direction of the bond.
    (See https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01424 Figure. 5)

    Args:
        bond_dipole: (unit: a.u. = e*A0)
            the norm of the bond dipole that use +/- to indicate
            the direction.
        bond_field_strength: (unit: kcal/mol*e*A)
            the EF projected on the bond. Need to use the
            same definition of the direction as bond_dipole
    Returns:
        g_ele = - E * u (unit: kcal/mol)
            the electrostatic stablization energy"""
    A0 = physical_constants["atomic unit of length"][0] # m
    A = 10**-10 # m
    g_ele = -bond_dipole * bond_field_strength * A0/A

    return g_ele

def ele_stab_energy_of_dipole(
        dipole: ArrayLike, 
        field_strength: ArrayLike,) -> float:
    """calculate the electrostatic stablization energy of a dipole.
    Not necessary for a bond.
    (See https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01424 Figure. 5)

    Args:
        dipole: (unit: a.u. = e*A0)
            the bond dipole as a vector.
        field_strength: (unit: kcal/mol*e*A)
            the EF as a vector in the same cartesian space as dipole

    Returns:
        g_ele = - E * u (unit: kcal/mol)
            the electrostatic stablization energy"""
    A0 = physical_constants["atomic unit of length"][0] # m
    A = 10**-10 # m
    g_ele = - np.dot(dipole, field_strength) * A0/A

    return g_ele

    
