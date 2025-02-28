"""Testing enzy_hpt.core.math_helper.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-26
"""
import logging
import numpy as np
import enzy_htp
from enzy_htp.core.general import EnablePropagate
from enzy_htp import _LOGGER
from enzy_htp.core import math_helper as mh
from scipy.spatial.transform import Rotation as R


def test_check_valid_ph_good_input(caplog):
    """Testing that the check_valid_ph() function works for good input."""
    existing_level = enzy_htp._LOGGER.level
    enzy_htp._LOGGER.setLevel(logging.DEBUG)
    phValues = np.linspace(0, 14, 100)
    with EnablePropagate(_LOGGER):
        for ph in phValues:
            mh.check_valid_ph(ph)
            assert "assigned pH:" not in caplog.text

    enzy_htp._LOGGER.setLevel(existing_level)


def test_check_valid_ph_bad_input(caplog):
    """Testing that the check_valid_ph() functio fails for bad input."""
    existing_level = _LOGGER.level
    _LOGGER.setLevel(logging.DEBUG)
    with EnablePropagate(_LOGGER):
        mh.check_valid_ph(-1)
        assert "assigned pH: -1.00 out of range: [0.00,14.00]" in caplog.text
        mh.check_valid_ph(15)
        assert "assigned pH: 15.00 out of range: [0.00,14.00]" in caplog.text
    enzy_htp._LOGGER.setLevel(existing_level)


def test_get_geom_center():
    """test function works as expected"""
    test_list_of_p = [(1, 1, 1), (2, 2, 5)]
    assert mh.get_geom_center(test_list_of_p) == (1.5, 1.5, 3)


def test_get_dihedral():
    """test if function works as expected"""
    test_points = [
        (33.49599838256836, 35.71099853515625, 26.554000854492188),
        (33.0880012512207, 37.07600021362305, 26.940000534057617),
        (34.00899887084961, 37.507999420166016, 28.075000762939453),
        (33.8380012512207, 38.87300109863281, 28.569000244140625),
    ]
    assert mh.get_dihedral(*test_points) == 176.4826686167864
    test_points = [
        (33.0880012512207, 37.07600021362305, 26.940000534057617),
        (34.00899887084961, 37.507999420166016, 28.075000762939453),
        (33.8380012512207, 38.87300109863281, 28.569000244140625),
        (34.38399887084961, 39.939998626708984, 28.013999938964844),
    ]
    assert mh.get_dihedral(*test_points) == -83.3566852742675
    test_points = [
        (39.07899856567383, 39.97600173950195, 37.40700149536133),
        (37.7400016784668, 39.50699996948242, 38.034000396728516),
        (36.72200012207031, 40.66899871826172, 37.95899963378906),
        (36.42900085449219, 41.4739990234375, 39.000999450683594),
    ]
    assert np.isclose(mh.get_dihedral(*test_points), -100.3235330367286, atol=1e-6)

def test_convert_first_three_point():
    """test using an example data"""
    test_internal_coord = [
        [0, 0.000, -1,    .0  , -2,     .0  ],
        [0, 1.449,  0,    .0  , -1,     .0  ],
        [1, 1.523,  0, 111.21 ,  0,     .0  ],
        [2, 1.540,  1, 111.208,  0, -180.000],
        [3, 1.218,  2, 107.632,  1, -179.951],
        [4, 1.218,  3, 119.996,  2,    0.229],
    ]
    result = mh._convert_first_three_point(test_internal_coord[:3])
    answer = [  [0., 0., 0.],
                [1.449, 0., 0.],
                [2.00000204, 1.419833, 0.],]
    for i,j in zip(result, answer):
        for k,l in zip(i,j):
            assert np.isclose(k,l)

def test_calcuate_cartesian():
    """test using an example data"""
    test_known_points = [
        np.array([0., 0., 0.]),
        np.array([1.449, 0., 0.]),
        np.array([2.00000204, 1.419833, 0.]),
    ]
    test_point = [2, 1.540,  1, 111.208,  0, 20.000]

    new_point =  mh._calculate_cartesian(test_point, test_known_points, False)
    answer = [0.94383, 2.42729, 0.49104]
    for i,j in zip(new_point, answer):
        assert np.isclose(i,j)

def test_internal_to_cartesian():
    """test using an example data.
    answer provided by http://www.shodor.org/chemviz/zmatrices/babel.html"""
    test_internal_coord = [
        [0, 0.000,  0,    .0  ,  0,     .0  ],
        [0, 1.449,  0,    .0  ,  0,     .0  ],
        [1, 1.523,  0, 111.21 ,  0,     .0  ],
        [2, 1.540,  1, 111.208,  0, -180.000],
        [3, 1.218,  2, 107.632,  1, -179.951],
        [4, 1.218,  3, 119.996,  2,    0.229],
    ]
    result = mh.internal_to_cartesian(test_internal_coord)
    answer = [
        [3.54000, 1.41978, -0.00000],
        [3.90898, 2.58055,  0.00079],
        [3.08816, 3.48041,  0.00562],
    ]
    for i,j in zip(result, answer):
        for k,l in zip(i,j):
            assert np.isclose(k,l, atol=1e-3)

def test_internal_to_cartesian_from_1():
    """test using an example data.
    answer provided by http://www.shodor.org/chemviz/zmatrices/babel.html"""
    test_internal_coord = [
        [0, 0.000,  0,    .0  ,  0,     .0  ],
        [1, 1.449,  0,    .0  ,  0,     .0  ],
        [2, 1.523,  1, 111.21 ,  0,     .0  ],
        [3, 1.540,  2, 111.208,  1, -180.000],
        [4, 1.218,  3, 107.632,  2, -179.951],
        [5, 1.218,  4, 119.996,  3,    0.229],
    ]
    result = mh.internal_to_cartesian(test_internal_coord, start_from_1=True)
    answer = [
        [3.54000, 1.41978, -0.00000],
        [3.90898, 2.58055,  0.00079],
        [3.08816, 3.48041,  0.00562],
    ]
    for i,j in zip(result, answer):
        for k,l in zip(i,j):
            assert np.isclose(k,l, atol=1e-3)

def test_get_angle():
    """test using an example data.
    answer get from pymol"""
    test_points = [
        (-46.808998107910156, 52.667999267578125, -28.253999710083008),
        (-47.63199996948242, 50.47800064086914, -28.195999145507812),
        (-47.42100143432617, 51.42599868774414, -30.475000381469727),
    ]
    assert np.isclose(mh.get_angle(*test_points), 65.7414, atol=1e-3)
    test_points = [
        (-20.92099952697754, 38.88800048828125, -11.484000205993652),
        (-20.18000030517578, 34.959999084472656, -12.663999557495117),
        (-22.488000869750977, 34.05500030517578, -12.204000473022461),
    ]
    assert np.isclose(mh.get_angle(*test_points), 97.116, atol=1e-3)

    test_points = [
        (0, 2, 0),
        (0, 0, 0),
        (0, 1, 0),
    ]
    assert np.isclose(mh.get_angle(*test_points), 0.0, atol=1e-3)

def test_is_integer():
    """as name"""
    test_num = 6.000001
    test_num_1 = 6.3
    test_num_2 = 5.99999
    assert mh.is_integer(test_num)
    assert not mh.is_integer(test_num_1)
    assert mh.is_integer(test_num_2)

def test_get_section_from_endpoint():
    """as name"""
    test_list = [(1,1),(10,2),(15,5),(25,5)]
    answer_list = [((1,1),(10,2)), ((11,2),(15,5)), ((16,5),(25,5))]
    assert mh.get_section_from_endpoint(test_list) == answer_list


def test_round_by():
    """as name"""

    assert mh.round_by(1.4, 0.5) == 1.0
    assert mh.round_by(-1.4, 0.5) == -1.0
    assert mh.round_by(-1.4, 0.3) == -2.0
    assert mh.round_by(1.4, 0.3) == 2.0


def test_rot_vec_from_dihedral():
    test_pt1 = np.array([1, 0, 0])
    test_pt2 = np.array([5, 7, 3])
    test_pt3 = np.array([1, 0, 5])
    test_pt4 = np.array([4, -5, 7])

    rv = mh.rot_vec_from_dihedral(test_pt1, test_pt2, test_pt3, test_pt4, 0)
    rot = R.from_rotvec(rv, degrees=True)

    dv = test_pt4 - test_pt3
 
    test_pt4 = test_pt3 + rot.apply(dv)

    assert abs(mh.get_dihedral(test_pt1, test_pt2, test_pt3, test_pt4)) <= 0.05

    test_pt1_2 = np.array([1, 0, 0])
    test_pt2_2 = np.array([5, 7, 3])
    test_pt3_2 = np.array([1, 0, 5])
    test_pt4_2 = np.array([4, -5, 7])

    rv2 = mh.rot_vec_from_dihedral(test_pt1_2, test_pt2_2, test_pt3_2, test_pt4_2, 120)
    rot2 = R.from_rotvec(rv2, degrees=True)

    dv2 = test_pt4_2 - test_pt3_2
 
    test_pt4_2 = test_pt3_2 + rot2.apply(dv2)

    assert abs(abs(mh.get_dihedral(test_pt1_2, test_pt2_2, test_pt3_2, test_pt4_2)) - 120.0) <= 0.05