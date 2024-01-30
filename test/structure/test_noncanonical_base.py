"""Testing the NonCanonicalBase() class in enzy_htp.structure. 
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-10-18
"""

from enzy_htp.structure import (
    NonCanonicalBase,
    Atom,
    )

def test_clone_connectivity():
    """test the function using fake objects"""
    # 123 are connected
    a1 = Atom.from_biopandas({"x_coord": 0, "y_coord": 0, "z_coord": 0, "atom_name": "O"})
    a2 = Atom.from_biopandas({"x_coord": 0.96, "y_coord": 0, "z_coord": 0, "atom_name": "D1"})
    a3 = Atom.from_biopandas({"x_coord": -0.24, "y_coord": 0.93, "z_coord": 1, "atom_name": "D2"})
    a1.connect = [(a2, "s"), (a3, "s")]
    a2.connect = [(a1, "s")]
    a3.connect = [(a1, "s")]
    # 456 are not connected
    a4 = Atom.from_biopandas({"x_coord": 0, "y_coord": 0, "z_coord": 0, "atom_name": "O"})
    a5 = Atom.from_biopandas({"x_coord": 0.96, "y_coord": 0, "z_coord": 0, "atom_name": "D1"})
    a6 = Atom.from_biopandas({"x_coord": -0.24, "y_coord": 0.93, "z_coord": 1, "atom_name": "D2"})

    ncaa_1 = NonCanonicalBase(
        residue_idx=1,
        residue_name="D2O",
        atoms=[a1, a2, a3],
    )
    ncaa_2 = NonCanonicalBase(
        residue_idx=1,
        residue_name="D2O",
        atoms=[a4, a5, a6],
    )

    ncaa_2.clone_connectivity(ncaa_1)

    assert ncaa_2.is_connected()
    assert ncaa_2.atoms[0].connect == [(a5, "s"), (a6, "s")]
    assert ncaa_2.atoms[1].connect == [(a4, "s")]
    assert ncaa_2.atoms[2].connect == [(a4, "s")]
