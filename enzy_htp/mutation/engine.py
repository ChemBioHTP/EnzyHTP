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
from .mutation import generate_all_mutations


def mutate_pdb(
    pdb: str,
    n_mutations: int = 1,
    mutations: List[str] = None,
    restrictions: MutationRestrictions = None,
    engine: str = "tleap",
    out_dir: str = None,
    random_state: int = 100,
) -> str:
    """Top-level method for mutating the structure contained in a .pdb file.

    Args:
        pdb: str() of the path to a .pdb file.
        n_mutations: number of desired mutations, default is 1.
        mutations: list() of supplied mutations.
        restrictions: MutationRestrictions() object which is created by calling TODO(CJ) on the input .pdb file.
        engine: The engine used to apply the mutation, default is 'tleap'.
        out_dir: The directory to save the mutated .pdb to. By default saves to the same directory as the base file.
        random_state: random state used for seeding mutation selection, default is 100.


    Returns:
        Path to the mutated .pdb file.

    Raises:
    """
    mutations = list(set(mutations))
    structure: struct.Structure = struct.structure_from_pdb(pdb)
    # TODO(CJ): need to update mutations
    # TODO(CJ): add warning if requested mutations are greater than
    # number of residues.
    if len(mutations) > n_mutations:
        _LOGGER.warn("TODO(CJ)")

    needed: int = n_mutations - len(mutations)
    excess: Dict[Tuple[str, int], List[Mutation]] = generate_all_mutations(structure)
    # TODO(CJ): remove the stuff we dont care about
    mutations.extend(get_n_mutations(excess, needed, restrictions, random_state))

    IMPLEMENTATION: Dict = {"tleap": apply_mutation_tleap}
    if engine not in IMPLEMENTATION:
        pass  # TODO(CJ): add some kind of error here
    return IMPLEMENTATION[engine](structure, mutations)


def get_n_mutations(
    excess: Dict[Tuple[str, int], List[Mutation]],
    needed: int,
    restrict: Union[MutationRestrictions, None],
    random_state: int,
):
    """ """
    np.random.seed(random_state)
    pass


def apply_mutations(
    structure: struct.Structure, mutations: List[Mutation], engine
) -> str:
    pass


def apply_mutation_tleap():
    pass
