"""EnzyHTP's physics knowledge.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-05-12
"""
from typing import Union, Tuple, List
import numpy as np
from numpy.typing import ArrayLike

from enzy_htp.core.logger import _LOGGER

def electric_field_strength(
        p0: ArrayLike,
        c0: float,
        p1: ArrayLike, 
        d1: ArrayLike = None,
        unit: str = "kcal/(mol*e*Ang)",
    ) -> Union[float, ArrayLike]:
    """return field strength E of source charge *p0(c0)* at *p1* in direction of *d1* if given
    Args:
        p0:
            position of charge of field source (unit: Ang)
        c0:
            point charge in p0 (unit: e)
        p1:
            the query position to calculate the field strength (unit: Ang)
        d1: (optional)
            the vector that defines the direction of field strength projection
        unit:
            the unit of the result (default: kcal/(mol*e*Ang))

    Returns:
        the field strength at the point (along the direction if provided)

    Details:
        -- E = kq/r^2 --
        point charge:   c0 in p0 
        point:          p1
        direction:      d1"""

    # Unit
    A = 10**-8 #cm
    Na = 6.02*10**23
    kcal = 4184 #J 
    e = 1.602*10**-19 #C     
    unit_scale_mapper = {
        "kcal/(mol*e*Ang)" : 1.0,
        "MV/cm" : kcal/(A*e*Na*10**6),
    }
    """map units to its scaling factor from `kcal/(mol*e*Ang)`"""

    # san check
    if unit not in unit_scale_mapper:
        _LOGGER.error(f"unit {unit} is not supported. (supported list: {unit_scale_mapper.keys()})")
        raise ValueError

    # init variables
    k = 332.4  # kcal*Ang/(mol*e^2) = (10^10)*(4.184^-1)*(Na)*(10^-3)*(1.602*10^-19)^2 * 9.0*10^-9 N*m^2/C^2
    q = c0  # e
    p0 = np.array(p0)  # Ang
    p1 = np.array(p1)  # Ang

    # Get r
    r = p1 - p0
    r_m = np.linalg.norm(r)  # magnitude of r
    r_u = r / r_m  # direction of r

    # Get E at p1
    e_ele_d = (k * q / (r_m**2)) * r_u # kcal/(mol*e*Ang)
    
    if d1 is not None:
        d1 = d1 / np.linalg.norm(d1)  # Ang
        # Get E along d1
        e_ele_d = np.dot(e_ele_d, d1)

    e_ele_d = e_ele_d * unit_scale_mapper[unit]

    return e_ele_d
