"""Module holds functions for math calculations. 

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-26
"""
from typing import Tuple, Union
from unittest.mock import NonCallableMock
from .logger import _LOGGER
import numpy as np
import math


def check_valid_ph(ph: float) -> NonCallableMock:
    """Helper function that checks the pH is on the range: [0.00, 14.00]. Give a warning if not."""
    if ph < 0 or ph > 14:
        _LOGGER.warning(f"assigned pH: {ph:.2f} out of range: [0.00,14.00]")


def set_distance(p1: Union[tuple, list], p2: Union[tuple, list],
                 d: float) -> Tuple[float, float, float]:
    """
    determine the coordinate of a point p3 with p1, p2, and d.
    p3 will be in the projection line of p1->p2 and have a d distance with p1.
    Args:
        p1: determine origin of the distance projection line
        p2: determine the direction of the distance projection line
        d : the length of the distance projection line
    Returns:
        the coord of p3
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    #direction vector
    v1 = (p2 - p1) / np.linalg.norm(p1 - p2)
    p3 = p1 + v1 * d

    return tuple(p3)


def get_distance(p1: Union[tuple, list], p2: Union[tuple, list]) -> float:
    """get the distance between p1 and p2"""
    d1 = np.array(p2) - np.array(p1)
    D = np.linalg.norm(d1)
    return D


def get_field_strength_value(p0: Union[tuple, list],
                             c0: float,
                             p1: Union[tuple, list],
                             p2: Union[tuple, list] = None,
                             d1: Union[tuple, list] = None) -> float:
    """
    calculate field strength E of *p0(c0)* at *p1* in direction of *p2-p1* or *d1*
    E = kq/r^2 (Unit: kcal/mol)
    point charge:   c0 in p0 
    point:          p1
    direction:      p2-p1 or d1
    Args:
        p0: position of charge of field source
        c0: point charge in p0
        p1: the query position to calculate the field strength
        p2: the point that defines the direction of field strength projection by p2-p1
        d1: the vector that defines the direction of field strength projection
    Returns:
        E: the field strength
    """
    # Unit
    k = 332.4  # kcal*Ang/(mol*e^2) = (10^10)*(4.184^-1)*(Na)*(10^-3)*(1.602*10^-19)^2 * 9.0*10^-9 N*m^2/C^2
    q = c0  # e
    p0 = np.array(p0)  # Ang
    p1 = np.array(p1)  # Ang
    if d1 == None:
        d1 = np.array(p2) - p1
    else:
        d1 = np.array(d1)
    d1 = d1 / np.linalg.norm(d1)  # Ang

    # Get r
    r = p1 - p0
    r_m = np.linalg.norm(r)
    r_u = r / r_m

    # Get E
    E = (k * q / (r_m**2)) * r_u
    # Get E in direction
    Ed = np.dot(E, d1)  # kcal/(mol*e*Ang)

    return Ed


def get_center(p1: Union[tuple, list],
               p2: Union[tuple, list]) -> Tuple[float, float, float]:
    """
    return the center of p1 and p2
    """
    p3 = 0.5 * (np.array(p1) + np.array(p2))

    return tuple(p3)


def round_by(num: float, cutnum: float) -> int:
    """
    round the float number up if the decimal part is larger than cutnum
    otherwise round down
    """
    dec_part, int_part = math.modf(num)
    if dec_part > cutnum:
        int_part += 1
    return int(int_part)
