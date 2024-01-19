"""Testing the functionality implemented in enzy_htp.chemical.physics.py

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-01-18"""
import numpy as np

import enzy_htp.chemical.physics as phy

def test_electric_field_strength():
    """as name. use result from internet as answer"""
    p0 = np.array((1.0, 1.0, 5.0))
    c0 = 2 #e
    p1 = np.array((10.0, 5.0, 2.0))
    d1 = np.array((3.0, 5.0, 7.0))
    answer_magn = 27.1691424 # MV/cm
    answer_direction_magn = np.dot((p1 - p0)/np.linalg.norm(p1 - p0), d1/np.linalg.norm(d1)) * answer_magn

    result = phy.electric_field_strength(
        p0, c0, p1, d1, unit = "MV/cm"
    )
    assert np.isclose(result, answer_direction_magn, atol=0.05)
