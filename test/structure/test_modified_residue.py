"""Testing the ModifiedResidue() class in enzy_htp.structure. 
Author: Sebastian Stull (sebastian.l.stull@vanderbilt.edu)
Date: 2025-02-18
"""

import os
from enzy_htp.structure import (
    ModifiedResidue,
    Atom,
    )
from enzy_htp.structure.structure_io.pdb_io import PDBParser

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURR_DIR}/data/"

def test_find_mainchain():

    a1 = Atom.from_biopandas({"x_coord": 0, "y_coord": 0, "z_coord": 0, "atom_name": "N"})
    a2 = Atom.from_biopandas({"x_coord": 0.96, "y_coord": 0, "z_coord": 0, "atom_name": "CA"})
    a3 = Atom.from_biopandas({"x_coord": -0.24, "y_coord": 0.93, "z_coord": 1, "atom_name": "C"})
    a4 = Atom.from_biopandas({"x_coord": -0.26, "y_coord": 1.94, "z_coord": 1, "atom_name": "CB"})
    a5 = Atom.from_biopandas({"x_coord": 0.26, "y_coord": -1.94, "z_coord": -1, "atom_name": "O"})

    # create fake connectivity
    a1.connect = [(a2, "s")]
    a2.connect = [(a1, "s"), (a3, "s")]
    a3.connect = [(a4, "s")]
    a4.connect = [(a5, "s")]
    a5.connect = [(a4, "s")]

    ncaa_1 = ModifiedResidue(
        residue_idx=1,
        residue_name="TST",
        atoms=[a1, a2, a3, a4, a5],
    )

    # assert main chain is correct
    mc = ncaa_1.find_mainchain()
    res = [aa.name for aa in mc]

    assert res == ["N", "CA", "C"]