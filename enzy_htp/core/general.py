"""This module contains helper functions related to general python/python library features
that does not suit any existing modules under the .core 
These functions will be grouped into specific modules in the future if they feel really
aggragates.

Current functions/classes in the file are:
List related:
    Class:
        GhostListElement
    Constant:
        GHOST_LIST_ELEMENT
    Functions:
        delete_base_on_id
        get_interval_from_list
        get_interval_str_from_list
        get_random_list_elem
        pop_random_list_elem
        product_lists_allow_empty

Dict related:
    Functions:
        get_copy_of_deleted_dict

Misc:
    Decorator:
        timer

Author: QZ Shao, <shaoqz@icloud.com>
Date: 2022-10-21
"""
import copy
import numpy as np
from typing import List, Iterable, Tuple, Dict
from time import perf_counter
import itertools

from .logger import _LOGGER


# == List related ==
class GhostListElement:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = object.__new__(cls)
            cls._instance._value = cls._instance
        return cls._instance

    def __eq__(self, other):
        return self is other


GHOST_LIST_ELEMENT = GhostListElement()
"""singleton for a ghost list element used for product_lists_allow_empty"""


def delete_base_on_id(target_list: list, target_id: int):
    """
    delete an element from a list base on its id() value
    """
    for i in range(len(target_list) - 1, -1, -1):
        if id(target_list[i]) == target_id:
            del target_list[i]


def get_interval_from_list(target_list: List[int]) -> Iterable[Tuple[int, int]]:
    """
    convert a list of int to the interval/range representation
    Returns:
        a generater of tuples with each indicating the start/end of the interval
    Example:
        >>> list(get_interval_from_list([1,2,3,6,7,8]))
        [(1,3),(6,8)]
    reference: https://stackoverflow.com/questions/4628333
    """
    # clean input
    target_list = sorted(set(target_list))
    # here use enum id as a ref sequence and group by the deviation
    for i, j in itertools.groupby(
            enumerate(target_list),
            lambda ref_vs_target: ref_vs_target[1] - ref_vs_target[0]):
        j = list(j)
        yield j[0][1], j[-1][1]


def get_interval_str_from_list(target_list: List[int]) -> str:
    """convert the result of get_interval_from_list() to a str"""
    result = get_interval_from_list(target_list)
    result = ",".join([f"{x[0]}-{x[1]}" if x[0] != x[1] else f"{x[0]}" for x in result])
    return result


def get_random_list_elem(traget_list: list):
    """Helper method that randomly chooses an element from a list. numpy.random.choice() doesn't 
    like to take elements from list()'s of tuples so this is the work around."""
    return np.random.choice(traget_list)


def pop_random_list_elem(traget_list: list):
    """Helper method that randomly pop an element from a list. (delete from original list)
    numpy.random.choice() doesn't like to take elements from list()'s of tuples so
    this is the work around."""
    return traget_list.pop(np.random.randint(len(traget_list)))


def product_lists_allow_empty(list_of_lists: List[list]) -> List[list]:
    """product lists in a list and give all possible cases of one element from each.
    Different from itertool.product, not getting any element from a list is allowed.
    Note the result will include an empty list"""
    list_of_lists_copy = copy.deepcopy(list_of_lists)
    list_of_lists_copy = [
        each_list + [GHOST_LIST_ELEMENT] for each_list in list_of_lists_copy
    ]
    return _product_lists_w_each_empty_ele(iter(list_of_lists_copy))


def _product_lists_w_each_empty_ele(list_of_lists: Iterable[list]) -> List[list]:
    """a sub-function used for product_list_allow_empty"""
    curr_list = next(list_of_lists, None)
    if not curr_list:
        return [[]]
    next_list = _product_lists_w_each_empty_ele(list_of_lists)
    return [[x] + y if x != GHOST_LIST_ELEMENT else y for x in curr_list
            for y in next_list]


# == Dict related ==
def get_copy_of_deleted_dict(orig_dict: Dict, del_key) -> Dict:
    """
    get a copy of the orig_dict and delete an item base on the del_key
    the del_key can be a list of keys
    """
    dict_copy = copy.deepcopy(orig_dict)
    if isinstance(del_key, list):
        for single_key in del_key:
            del dict_copy[single_key]
    else:
        del dict_copy[del_key]

    return dict_copy


# == misc ===


def timer(fn):
    """decodator for timing the run of the function {fn}"""

    def inner(*args, **kwargs):
        start_time = perf_counter()
        to_execute = fn(*args, **kwargs)
        end_time = perf_counter()
        execution_time = end_time - start_time
        _LOGGER.info("{0} took {1:.8f}s to execute".format(fn.__name__, execution_time))
        return to_execute

    return inner
