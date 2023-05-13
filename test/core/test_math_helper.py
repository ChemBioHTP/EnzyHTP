"""Testing enzy_hpt.core.math_helper.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-26
"""
import logging
import numpy as np
import enzy_htp
from enzy_htp.core import math_helper as mh


def test_check_valid_ph_good_input(caplog):
    """Testing that the check_valid_ph() function works for good input."""
    existing_level = enzy_htp._LOGGER.level
    enzy_htp._LOGGER.setLevel(logging.DEBUG)
    phValues = np.linspace(0, 14, 100)
    for ph in phValues:
        mh.check_valid_ph(ph)
        assert "assigned pH:" not in caplog.text

    enzy_htp._LOGGER.setLevel(existing_level)


def test_check_valid_ph_bad_input(caplog):
    """Testing that the check_valid_ph() functio fails for bad input."""
    existing_level = enzy_htp._LOGGER.level
    enzy_htp._LOGGER.setLevel(logging.DEBUG)
    mh.check_valid_ph(-1)
    assert "assigned pH: -1.00 out of range: [0.00,14.00]" in caplog.text
    mh.check_valid_ph(15)
    assert "assigned pH: 15.00 out of range: [0.00,14.00]" in caplog.text
    enzy_htp._LOGGER.setLevel(existing_level)

def test_get_geom_center():
    """test function works as expected"""
    test_list_of_p = [
        (1, 1, 1),
        (2, 2, 5)
    ]
    assert mh.get_geom_center(test_list_of_p) == (1.5, 1.5, 3)
