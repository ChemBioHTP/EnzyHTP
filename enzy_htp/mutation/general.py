"""Define functions for mutate Structure(). 
Science API:
    +mutate_stru()

Mutation is carried out by an underlying engine and the supported engines currently include:
    + Amber/tleap

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-10-24
"""

from enzy_htp.structure import Structure

def mutate_stru(stru: Structure,
                pattern: str, # because we need to support insertion or deletion
                engine: str = "tleap",
                in_place: bool=False) -> Structure:
    """
    This science API solves the protein mutation problem.
    It means to determine the structure of the mutant of the source protein caused
    by residue substitution, deletion, and insertion. (see also: structure prediction
    for alternative functions)
    Args:
        stru: the target structure
        pattern: the pattern for describing a set of mutations to deploy
        in_place: if change the structure in-place and return the reference
                  default is False since wild-type structure is expected to also available
                  in many applications
    Returns:
        the reference/copy of the changed structure
    """

# mutate stru
# assign mutation
# deploy mutation