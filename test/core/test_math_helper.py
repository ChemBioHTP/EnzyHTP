"""Testing enzy_hpt.core.math_helper.py

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-26
"""
import numpy as np
from enzy_htp.core import math_helper as mh

def test_check_valid_ph_good_input(caplog):
    """Testing that the check_valid_ph() functio works for good input."""
    phValues = np.linspace(0, 14, 100)
    for ph in phValues:
        mh.check_valid_ph(ph)
        assert "assigned pH:" not in caplog.text

def test_check_valid_ph_bad_input(caplog):
    """Testing that the check_valid_ph() functio fails for bad input."""
    mh.check_valid_ph(-1)
    assert "assigned pH: -1.00 out of range: [0.00,14.00]" in caplog.text
    mh.check_valid_ph(15)
    assert "assigned pH: 15.00 out of range: [0.00,14.00]" in caplog.text