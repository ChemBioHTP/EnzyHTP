"""EnzyHTP's physics knowledge.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-05-12
"""
from typing import Union, Tuple, List
import numpy as np

def get_ele_field_strength_value(p0: Union[Tuple, List],
                                 c0: float,
                                 p1: Union[Tuple, List],
                                 p2: Union[Tuple, List] = None,
                                 d1: Union[Tuple, List] = None) -> float:
    '''
    return field strength E of *p0(c0)* at *p1* in direction of *p2-p1* or *d1*
    Details:
        -- E = kq/r^2 -- (Unit: kcal/(mol*e*Ang))
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
        e_ele_d: the field strength on the direction
    '''
    # Unit
    k = 332.4   # kcal*Ang/(mol*e^2) = (10^10)*(4.184^-1)*(Na)*(10^-3)*(1.602*10^-19)^2 * 9.0*10^-9 N*m^2/C^2
    q = c0                      # e
    p0 = np.array(p0)           # Ang
    p1 = np.array(p1)           # Ang
    if d1 is None:
        d1 = np.array(p2) - p1
    else:
        d1 = np.array(d1)
    d1 = d1/np.linalg.norm(d1)  # Ang

    # Get r
    r = p1 - p0
    r_m = np.linalg.norm(r)
    r_u = r/r_m

    # Get E
    e_ele = (k * q / (r_m**2)) * r_u
    # Get E in direction
    e_ele_d = np.dot(e_ele, d1)          # kcal/(mol*e*Ang)

    return e_ele_d

