"""Sub-module in enzy_htp.chemical holding information about chemical bond lengths.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-16
"""
from .db import load_from_db

from typing import Dict

X_H_BOND_LENGTH: Dict[str, float] = load_from_db('X_H_BOND_LENGTH')
""""""


def get_h_bond_length(aname: str) -> float:
    """"""
    pass
