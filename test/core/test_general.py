"""Testing enzy_hpt.core.general.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-02-03
"""

import os
import logging
import numpy as np
from enzy_htp.core import general as eg
from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data"
WORK_DIR = f"{CURR_DIR}/work_dir"

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

def test_split_but_brackets():
    """as name"""
    test_pattern = "'B.254.CAE', 'B.254.H2', (0,10)"
    assert len(eg.split_but_brackets(test_pattern, ",")) == 3

    test_pattern = "[1,1], {2,2}, (0,10)"
    assert len(eg.split_but_brackets(test_pattern, ",")) == 3

    test_pattern = "[(1,1)], (0,10)"
    assert len(eg.split_but_brackets(test_pattern, ",")) == 2

    test_pattern = "'[1,1, [1,1]], (0,10,(1,1))"
    assert len(eg.split_but_brackets(test_pattern, ",")) == 2

def test_split_but_brackets_more_layer():
    """as name. a failed case. TODO"""
    test_pattern = "'[1,1, [1,1], 1], (0,10,(1,1))"
    assert len(eg.split_but_brackets(test_pattern, ",")) == 2

def test_load_obj():
    """test load_obj()"""
    test_pickle_file = f"{DATA_DIR}/test_multi.pickle"
    result = eg.load_obj(test_pickle_file)
    assert isinstance(result, list)
    assert len(result) == 3

    test_pickle_file = f"{DATA_DIR}/test_single.pickle"
    result = eg.load_obj(test_pickle_file)
    assert isinstance(result, int)

def test_save_obj():
    """test save_obj()"""
    test_pickle_file = f"{WORK_DIR}/test_save.pickle"
    test_obj = range(10)
    eg.save_obj(test_obj, test_pickle_file)
    answer_obj = eg.load_obj(test_pickle_file)
    assert answer_obj == test_obj

    # in case the file exists
    test_obj_1 = range(20)
    eg.save_obj(test_obj_1, test_pickle_file)
    answer_obj = eg.load_obj(test_pickle_file)
    assert answer_obj == test_obj_1

    fs.clean_temp_file_n_dir([test_pickle_file])


def test_log_level():
    """test LogLevel context manager"""
    # Create a test logger
    test_logger = logging.getLogger("test_logger")
    original_level = logging.INFO
    test_logger.setLevel(original_level)
    
    # Test that the context manager changes and restores the log level
    with eg.LogLevel(test_logger, logging.DEBUG):
        assert test_logger.level == logging.DEBUG
    
    # Check that the original level is restored
    assert test_logger.level == original_level
    
    # Test with different level
    with eg.LogLevel(test_logger, logging.ERROR):
        assert test_logger.level == logging.ERROR
    
    # Check that the original level is restored again
    assert test_logger.level == original_level


def test_log_level_with_exception():
    """test LogLevel context manager when exception occurs"""
    test_logger = logging.getLogger("test_logger_exception")
    original_level = logging.WARNING
    test_logger.setLevel(original_level)
    
    # Test that the original level is restored even when exception occurs
    try:
        with eg.LogLevel(test_logger, logging.DEBUG):
            assert test_logger.level == logging.DEBUG
            raise ValueError("Test exception")
    except ValueError:
        pass
    
    # Check that the original level is restored after exception
    assert test_logger.level == original_level
