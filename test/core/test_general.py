"""Testing enzy_hpt.core.general.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-02-03
"""

import numpy as np
from enzy_htp.core import general as eg
from enzy_htp.core import _LOGGER

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


def test_list_remove_adjacent_duplicates():
    """test function works as expected"""
    test_list = [1, 2, 2, 3, 3, 4, 6, 6, 3, 5, 5]
    result_list = eg.list_remove_adjacent_duplicates(test_list)
    assert result_list == [1, 2, 3, 4, 6, 3, 5]

def test_swapped_dict():
    """test function works as expected"""
    test_dict = {
        1 : "1",
        2 : "2",
        3 : "3",
    }
    result_dict = eg.swapped_dict(test_dict)
    assert result_dict == {
        "1" : 1,
        "2" : 2,
        "3" : 3,
    }

def test_num_ele_2d():
    """as name"""
    test_list_2d = [
        [1,2,3],
        [4,5],
        [6,7,8,9]
    ]
    assert eg.num_ele_2d(test_list_2d) == 9

def test_capture_logging(capfd):
    """as name"""
    with eg.CaptureLogging(_LOGGER) as log_str:
        _LOGGER.error("redirect test")
    
    assert log_str.getvalue() == "redirect test\n"
    # test handlers restored
    _LOGGER.error("restore test")

    captured = capfd.readouterr()
    assert "redirect" not in captured.err
    assert "restore" in captured.err
