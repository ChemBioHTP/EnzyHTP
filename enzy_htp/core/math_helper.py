"""Module holds functions for math calculations. 

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-26
"""
from typing import List, Tuple, Union
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

def get_dihedral(p1: Union[tuple, list], p2: Union[tuple, list],
                 p3: Union[tuple, list], p4: Union[tuple, list],
                 determine_sign: bool= True,
                 rad_result: bool= False) -> float:
    """get the dihedral defined by p1, p2, o3, and p4
    Args:
        p1, p2, p3, p4:
            the coordinate of the 4 defining points in tuple or list
            the 2 plane is defined by p1, p2, p3 and p2, p3, p4.
        determine_sign: (default: True)
            whether want to determine the sign of the result, if do:
                plane(123) spin to plane(234), watching along z-->z+,
                clockwise: +
                counter clockwise: -
        rad_result: (default: False)
            whether giving result as radian
    Return:
        the dihedral value (degrees by defaul)"""
    v1 = np.array(p2) - np.array(p1)
    v2 = np.array(p3) - np.array(p2)
    v3 = np.array(p4) - np.array(p3)
    nv1 = np.cross(v1, v2)
    nv2 = np.cross(v2, v3)
    dot = np.dot(nv1, nv2)
    mag1 = np.linalg.norm(nv1)
    mag2 = np.linalg.norm(nv2)
    dihedral = np.arccos(dot / (mag1 * mag2))

    if determine_sign:
        if np.cross(nv1, nv2)[2] < 0:
            dihedral *= -1
    if not rad_result:
        dihedral = np.degrees(dihedral)

    return dihedral


def get_center(p1: Union[tuple, list],
               p2: Union[tuple, list]) -> Tuple[float, float, float]:
    """
    return the center of p1 and p2
    """
    p3 = 0.5 * (np.array(p1) + np.array(p2))

    return tuple(p3)

def get_geom_center(list_of_p: list) -> Tuple[float, float, float]:
    """get the geometric center of a list of points"""
    return tuple(np.mean(np.array(list_of_p), axis=0))

def round_by(num: float, cutnum: float) -> int:
    """
    round the float number up if the decimal part is larger than cutnum
    otherwise round down
    """
    dec_part, int_part = math.modf(num)
    if dec_part > cutnum:
        int_part += 1
    return int(int_part)

def calc_average_task_num(num_of_task: int, num_of_worker: int) -> List[int]:
    """calculate task number for each worker based on {num_of_task} and
    {num_of_worker}. Return a list of task numbers for each worker."""
    task_each_worker = num_of_task // num_of_worker
    remaining_tasks = num_of_task % num_of_worker
    result = [task_each_worker for i in range(num_of_worker)]
    for i in range(remaining_tasks):
        result[i] += 1
    return result
