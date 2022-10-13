"""Testing the constants found in enzy_htp.chemical.metal

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-19
"""
import pytest
from typing import List

from enzy_htp.chemical import metal as mm

from util import all_caps

IONIC_ATOM_NAMES: List[str] = [
    "N",
    "O",
    "S",
    "Mg",
    "Li",
    "Zn",
    "Co",
    "Fe",
    "Mn",
    "Ca",
]

VDW_ATOM_NAMES: List[str] = [
    "N",
    "O",
    "S",
    "Mg",
    "Li",
    "Tl",
    "Cu",
    "Ag",
    "Be",
    "Ni",
    "Pt",
    "Zn",
    "Pd",
    "Hg",
    "Cd",
    "Ca",
    "Sn",
    "Pb",
    "Sr",
    "Ba",
    "Ra",
    "Al",
    "In",
]


def test_variable_checks():
    """Basic testing for the metal mapper found in enzy_htp.chemical.metal"""
    mmapper_keys = list(mm.METAL_MAPPER.keys())
    mmapper_keys.remove("Na+")
    assert all_caps(mmapper_keys)
    assert max(list(map(len, mmapper_keys)))

    mcenter_keys = list(mm.METAL_CENTER_MAP.keys())
    assert all_caps(mcenter_keys)
    assert max(list(map(len, mcenter_keys)))


def test_get_metal_radii_vdw():
    """Making sure the get_metal_radii() method works for all values in the original VDW_RADII dict()."""

    for aname in VDW_ATOM_NAMES:
        assert mm.get_atom_radii(aname, "vdw") is not None


def test_get_metal_radii_ionic():
    """Making sure the get_metal_radii() method works for all values in the original IONIC_RADII dict()."""

    for aname in IONIC_ATOM_NAMES:
        assert mm.get_atom_radii(aname, "ionic") is not None


def test_get_metal_radii_bad_inputs():
    """Making sure the get_metal_radii() method returns None for bad inputs both non-existent and case-insensitive."""

    for aname in VDW_ATOM_NAMES:
        if len(aname) < 2:
            continue

        with pytest.raises(SystemExit) as exe:
            mm.get_atom_radii(aname.upper(), "vdw")

        assert exe.type == SystemExit
        assert exe.value.code == 1

    for aname in IONIC_ATOM_NAMES:
        if len(aname) < 2:
            continue

        with pytest.raises(SystemExit) as exe:
            mm.get_atom_radii(aname.upper(), "ionic")

        assert exe.type == SystemExit
        assert exe.value.code == 1

    with pytest.raises(SystemExit) as exe:
        mm.get_atom_radii("DNE")

    assert exe.type == SystemExit
    assert exe.value.code == 1

    with pytest.raises(SystemExit) as exe:
        mm.get_atom_radii("C")

    assert exe.type == SystemExit
    assert exe.value.code == 1

    with pytest.raises(SystemExit) as exe:
        mm.get_atom_radii("H")

    assert exe.type == SystemExit
    assert exe.value.code == 1


def test_get_metal_radii_bad_input():
    """Ensuring that get_metal_radii() exits for an invalid input."""
    with pytest.raises(SystemExit) as exe:
        mm.get_atom_radii("Li", "not-real")

    assert exe.type == SystemExit
    assert exe.value.code == 1
