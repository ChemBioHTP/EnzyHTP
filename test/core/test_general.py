"""Testing enzy_hpt.core.general.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-02-03
"""

import numpy as np
from enzy_htp.core import general as eg


def test_pop_random_list_elem():
    """test function works as expected"""
    test_list = [1, 2, 3]
    np.random.seed(3)
    popped_elem = eg.pop_random_list_elem(test_list)
    assert popped_elem == 3
    assert test_list == [1, 2]


def test_product_lists_allow_empty():
    """test function works as expected"""
    test_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    result_list = eg.product_lists_allow_empty(test_list)
    assert len(result_list) == 64
    assert [1] in result_list
    assert [] in result_list
    assert test_list == [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
