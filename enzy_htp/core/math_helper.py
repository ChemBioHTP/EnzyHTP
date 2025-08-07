"""Module holds functions for math calculations. 

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-26
"""
from typing import List, Tuple, Union
from unittest.mock import NonCallableMock
from .logger import _LOGGER
import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation as R
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


def get_angle(p1: Union[tuple, list],
              p2: Union[tuple, list],
              p3: Union[tuple, list],
              rad_result: bool = False) -> float:
    """get the angle between p1, p2, p3.
    result range is [0,180]"""
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    v21 = p1 - p2
    v23 = p3 - p2
    result = get_angle_vec(v21, v23, rad_result)

    return result


def get_angle_vec(
        d1: Union[tuple, list],
        d2: Union[tuple, list],
        rad_result: bool = False) -> float:
    """get the angle between d1, d2
    result range is [0,180]"""
    n1 = np.cross(d2, d1)
    # consider collinear
    if np.isclose(np.linalg.norm(n1), 0, atol=1e-9):
        if np.dot(d1, d2) < 0:
            result = np.pi
        else:
            result = 0.0
    else:
        x_unit = d2 / np.linalg.norm(d2)
        y_unit = np.cross(n1, d2)
        y_unit = y_unit / np.linalg.norm(y_unit)
        x = np.dot(d1, x_unit)
        y = np.dot(d1, y_unit)

        result = np.arctan2(y, x)

    if not rad_result:
        result = np.degrees(result)

    return result


def get_dihedral(p1: Union[tuple, list],
                 p2: Union[tuple, list],
                 p3: Union[tuple, list],
                 p4: Union[tuple, list],
                 rad_result: bool = False) -> float:
    """get the dihedral defined by p1, p2, o3, and p4
    Args:
        p1, p2, p3, p4:
            the coordinate of the 4 defining points in tuple or list
            the 2 plane is defined by p1, p2, p3 and p2, p3, p4.
            plane(123) spin to plane(234), watching along z-->z+,
            clockwise: +
            counter clockwise: -
            ref: https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
        rad_result: (default: False)
            whether giving result as radian
    Return:
        the dihedral value (degrees by defaul)"""
    v1 = np.array(p2) - np.array(p1)
    v2 = np.array(p3) - np.array(p2)
    v3 = np.array(p4) - np.array(p3)
    nv1 = np.cross(v1, v2)
    nv1 = nv1 / np.linalg.norm(nv1)
    nv2 = np.cross(v2, v3)
    nv2 = nv2 / np.linalg.norm(nv2)

    m1 = np.cross(nv1, v2 / np.linalg.norm(v2))
    x = np.dot(nv1, nv2)
    y = np.dot(m1, nv2)
    dihedral = -np.arctan2(y, x)

    if not rad_result:
        dihedral = np.degrees(dihedral)

    return dihedral


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
    if abs(dec_part) > cutnum:
        if dec_part > 0:
            int_part += 1
        else:
            int_part -= 1
    return int(int_part)


def is_integer(num: float, tolerance: float = 0.01) -> bool:
    """determine whether a float is actually an integer within an error tolerance.
    The tolerence must be smaller than 1.
    Example: 1.999 is 2 by tolerance of 0.01"""
    # san check
    if tolerance > 1:
        _LOGGER.error(f"tolerance must be smaller than 1. (input: {tolerance})")
        raise ValueError

    dec_part, int_part = math.modf(num)
    return (dec_part < tolerance) or (dec_part > (1 - tolerance))


def get_section_from_endpoint(endpoint_list: List[Tuple[int, int]]) -> List[Tuple[Tuple[int, int], Tuple[int, int]]]:
    """get sections that does not share endpoints from a list of raw endpoints.
    e.g.: [(1,1),(10,2),(15,5)] -> [((1,1),(10,2)), ((11,2),(15,2))]"""
    result = [(endpoint_list[0], endpoint_list[1])]
    last_point = endpoint_list[1]
    for endpoint in endpoint_list[2:]:
        x, y = last_point
        non_share_point = (x+1, y)
        result.append((non_share_point, endpoint))
        last_point = endpoint
    return result


def calc_average_task_num(num_of_task: int, num_of_worker: int) -> List[int]:
    """calculate task number for each worker based on {num_of_task} and
    {num_of_worker}. Return a list of task numbers for each worker."""
    task_each_worker = num_of_task // num_of_worker
    remaining_tasks = num_of_task % num_of_worker
    result = [task_each_worker for i in range(num_of_worker)]
    for i in range(remaining_tasks):
        result[i] += 1
    return result


def internal_to_cartesian(internal_coordinate: List[npt.ArrayLike],
                          remove_dummy_point: bool=True,
                          start_from_1: bool=False,) -> List[npt.NDArray]:
    """derive cartesian coordinate for a list of points in internal coordinate.
    Args:
        internal_coordinate:
            format: [(bond_index, bond, angle_index, angle, dihedral_index, dihedral), ...]
        remove_dummy_point:
            whether remove the first 3 dummy points from the result
        start_from_1:
            whether the point indexes from {internal_coordinate} starts from 1
            instead of 0.
    Returns:
        cartesian_coordinate:
            format: [(x,  y, z), ...]"""
    result = []
    result.extend(_convert_first_three_point(internal_coordinate[:3]))

    for point in internal_coordinate[3:]:
        result.append(_calculate_cartesian(point, result, start_from_1))

    if remove_dummy_point:
        result = _remove_dummy_point(result)

    return result


def _convert_first_three_point(frist_three: List[npt.ArrayLike]) -> List[npt.NDArray]:
    """convert the first three point of internal coordinate using a different method"""
    result = []
    # 1st point
    p_1 = np.array([0.0, 0.0, 0.0])
    result.append(p_1)

    # 2nd point
    p_2 = np.array([frist_three[1][1], 0.0, 0.0])
    result.append(p_2)

    # 3rd point
    r = frist_three[2][1]
    theta = frist_three[2][3]
    p_3 = p_2 + np.array([r * np.cos(np.deg2rad(180.0-theta)), r * np.sin(np.deg2rad(180.0-theta)), 0])
    result.append(p_3)

    return result


def _calculate_cartesian(point: List,
                         known_points: List[npt.NDArray],
                         start_from_1: bool) -> npt.NDArray:
    """calculate the cartesian coordinate of {point} in internal coordinate
    using {known_points} in cartesian coordinate"""
    p1_id, r, p2_id, angle, p3_id, dihedral = point
    if start_from_1:
        p1_id -= 1
        p2_id -= 1
        p3_id -= 1

    p1 = known_points[p1_id]
    p2 = known_points[p2_id]
    p3 = known_points[p3_id]

    # calculate the rotation vector (unit_axis * theta_in_degree) for angle
    v1 = p2 - p3
    v2 = p1 - p2
    nv1 = np.cross(v1, v2)
    rv1 = nv1 / np.linalg.norm(nv1) * -angle
    rot_1 = R.from_rotvec(rv1, degrees=True)
    # calculate the rotation vector (unit_axis * theta_in_degree) for dihedral
    rv2 = v2 / np.linalg.norm(v2) * dihedral
    rot_2 = R.from_rotvec(rv2, degrees=True)

    # calculate displacement vector
    dv = r * (p2 - p1) / np.linalg.norm(p2 - p1)

    # rotate dv with rv1, rv2
    dv = rot_2.apply(rot_1.apply(dv))

    # apply dv to p1
    result = p1 + dv

    return result


def rotation_matrix_from_vectors(vec1, vec2): # TODO complete type hinting
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    
    NOTE: probably should use R.from_rotvec instead
    """ # TODO make unit test
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def _remove_dummy_point(coord: List[npt.NDArray]) -> List[npt.NDArray]:
    """remove the first 3 dummy points"""
    return coord[3:]


def direction_unit_vector(p1:npt.NDArray, p2:npt.NDArray) ->  npt.NDArray:
    """Calculates the unit vector from p1 to p2, where p1 and p2 are both.  The resulting
    vector "points" from p1 to p2. Requires that both p1 and p2 are numpy arrays."""

    result = p2 - p1

    return result / np.linalg.norm(result)

def rot_vec_from_dihedral(p0:npt.NDArray, p1:npt.NDArray, p2:npt.NDArray, p3:npt.NDArray, 
                          value:float) -> npt.NDArray:
    """Provides a rotation vector that sets the dihedral of the four given points 
    to a given value in the specified direction."""
    dihedral = get_dihedral(p0, p1, p2, p3)
    direction = direction_unit_vector(p1, p2)
    v = value - dihedral
    return (direction * v)