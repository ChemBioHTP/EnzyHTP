"""Sub-module in enzy_htp.chemical holding information about chemical bond lengths.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-16
"""
from .db import load_from_db

from typing import Dict

X_H_BOND_LENGTH: Dict[str, float] = load_from_db("X_H_BOND_LENGTH")
"""dict() storing bond lengths between an arbitrary atom X and a Hydrogen (H)."""


def get_h_bond_length(aname: str) -> float:
    """Gets the length of a bond with Hydrogen with an arbitrary atom, returning an answer in angstroms.

    Args:
        aname: The name of the atom to find the Hydrogen bond length for.

    Returns:
        The appropriate bond length in angstroms, or -1.0 if the value is not defined.
    """
    return X_H_BOND_LENGTH.get(aname, -1.0)
