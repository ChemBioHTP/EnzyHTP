"""Sub-module in enzy_htp.chemical holding information about chemical bond lengths.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-16
"""

from typing import Dict, List

X_H_BOND_LENGTH: Dict[str, float] = {
        "C": 1.07,
        "N": 1.0
}
"""dict() storing bond lengths between an arbitrary atom X and a Hydrogen (H)."""


def get_h_bond_length(aname: str) -> float:
    """Gets the length of a bond with Hydrogen with an arbitrary atom, returning an answer in angstroms.

    Args:
        aname: The name of the atom to find the Hydrogen bond length for.

    Returns:
        The appropriate bond length in angstroms, or -1.0 if the value is not defined.
    """
    return X_H_BOND_LENGTH.get(aname, -1.0)


def get_valid_generic_atom_name(atom_names: List[str]) -> List[str]:
    """
    Get valid generic atom names for a list of atom names. Normally this list of atom names
    come from a residue and the goal is to avoid repeating names. new names will add index from 1
    if any repeating ones are found.
    List index will hold the same to make the correpondance
    """
    result = []
    for name in atom_names:
        new_name = name
        count = 0
        while new_name in result:
            count += 1
            new_name = name + str(count)
        result.append(new_name)
    return result
