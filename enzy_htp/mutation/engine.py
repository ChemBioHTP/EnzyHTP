"""Top level engine for mutating the structure in a .pdb file.


Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-15
"""
import numpy as np
from typing import List, Dict


from enzy_htp.core import _LOGGER
import enzy_htp.structure as struct
from .mutation_restrictions import MutationRestrictions


def mutate_pdb(
    pdb: str,
    out_dir: str = None,
    n_mutations: int = 1,
    mutations: List[str] = None,
    restrictions: MutationRestrictions = None,
    engine: str = "tleap",
    random_state: int = 100,
) -> str:
    """Top-level method for mutating the structure contained in a .pdb file.

    Args:
        pdb: str() of the path to a .pdb file.
        out_dir: The directory to save the mutated .pdb to. By default saves to the same directory as the base file.
        n_mutations: number of desired mutations, default is 1.
        mutations: list() of supplied mutations.
        restrictions:
        engine: The engine used to apply the mutation, default is 'tleap'.
        random_state: random state used for seeding mutation selection, default is 100.


    Returns:
        Path to the mutated .pdb file.

    Raises:
    """
    mutations = list(set(mutations))
    structure: struct.Structure = struct.structure_from_pdb(pdb)
    if len(mutations) > n_mutations:
        _LOGGER.warn("TODO(CJ)")
    else:
        needed: int = n_mutations - len(mutations)
        mutations.extend(
            get_n_mutations(structure, needed, restrictions, mutations, random_state)
        )

    return apply_mutations(structure, mutations, engine)


def get_n_mutations(
    structure: struct.Structure,
    n: int,
    restrict: Dict,
    existing: List[str],
    random_state: int,
):
    """"""
    np.random.seed(random_state)
    pass


def apply_mutations(structure: struct.Structure, mutations, engine) -> str:
    pass


def valid_mutation(structure: struct.Structure, mutation: str, restriction) -> bool:
    """Checks if a given mutation is valid for the given Structure and set of restrictions.

    Args:
            structure: The Structure() object mutations will be performed on.
            mutation: The given mutation str() that will be applied to the Structure().
    """
    pass
