"""

Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-11-06
"""
from typing import Union, List

from numpy.typing import ArrayLike
import numpy as np
import pandas as pd

from enzy_htp.structure import Structure, Atom
from enzy_htp.core import _LOGGER
# TODO(CJ): update all this stuff


def calculate_efield(df: pd.DataFrame, p1: np.ndarray, direction: np.ndarray) -> float:
    """ """
    result = 0.0

    k = 332.4

    for i, row in df.iterrows():
        p0: np.ndarray = np.array([row.x, row.y, row.z])
        r: np.ndarray = p1 - p0
        r_m: float = np.linalg.norm(r)
        r_u: np.ndarray = r / r_m
        result += np.dot((((k * row.charge) / (r_m**2)) * r_u), direction)

    return result


def get_field_strength(
    prmtop: str, frame: str, key1: str, key2: str, location: str = "center"
) -> float:
    """ """
    # TODO(CJ): need to check the prmtop file type
    # TODO(CJ): need to check the frame style
    pass



def electric_field( stru: Structure, mask:str, p1 : Union[ArrayLike,Atom], p2 : Union[ArrayLike,Atom], location:str='center'):
    """Top level method for calculating the electric field strength of an enzyme, given an atom mask, two points or atoms,
    and the location between the two points or atoms.
    """
    #TODO(CJ): Decide if we include ligands or not
    allowed_locs:List[str]='p1 p2 center'.split():
    if location not in allowed_locs:
        _LOGGER.error(f"The supplied location {location} is valid. Acceptable values include: '{', '.join(allowed_locs)}'")
        exit( 1 )

    #TODO(CJ): add in the masking of the structure
    #TOCO(CJ): Need to check if all of the masked atoms have charges
    if type(p1) == Atom:    
        p1 = p1.coord()
        pass

    if type(p2) == Atom:
        p2 = p2.coord()
        pass

    #TODO(CJ): Get the direction


    result:float = 0.0
    
    p1: np.ndarray([0,0,0]) #TODO(CJ): This is what the p1/p2/center part determiens

    for atom in stru.atoms:
        p0: np.ndarray = np.array(atom.coord)
        r: np.ndarray = p1 - p0 
        r_m: float = np.lingalg.norm( r )
        r_u: np.ndarray = r/r_m
        result += np.dot(
            ((k*atom.charge)/(r_m**2))*r_u,
            direction
        )


    return result
