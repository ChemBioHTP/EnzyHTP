"""This module contains helper functions related to general python/python library features

Author: QZ Shao, <shaoqz@icloud.com>
Date: 2022-10-21
"""
import copy
import numpy as np
from typing import List, Iterable, Tuple, Dict
import itertools

# == List related ==
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
    return traget_list[np.random.randint(len(traget_list))]

def pop_random_list_elem(traget_list: list):
    """Helper method that randomly pop an element from a list. (delete from original list)
    numpy.random.choice() doesn't like to take elements from list()'s of tuples so
    this is the work around."""
    return traget_list.pop(np.random.randint(len(traget_list)))

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
