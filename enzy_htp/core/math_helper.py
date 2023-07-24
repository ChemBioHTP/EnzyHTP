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


def set_distance(p1: Union[tuple, list], p2: Union[tuple, list], d: float) -> Tuple[float, float, float]:
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


def get_center(p1: Union[tuple, list], p2: Union[tuple, list]) -> Tuple[float, float, float]:
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
