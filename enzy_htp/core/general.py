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
from io import StringIO
import os
import re
import sys
import logging
import time
import numpy as np
from typing import Any, List, Iterable, Tuple, Dict, Callable
import itertools
import pickle
import inspect

from .logger import _LOGGER
from .file_system import write_lines

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
    for i, j in itertools.groupby(enumerate(target_list), lambda ref_vs_target: ref_vs_target[1] - ref_vs_target[0]):
        j = list(j)
        yield j[0][1], j[-1][1]


def get_interval_str_from_list(target_list: List[int]) -> str:
    """convert the result of get_interval_from_list() to a str"""
    result = get_interval_from_list(target_list)
    result = ",".join([f"{x[0]}-{x[1]}" if x[0] != x[1] else f"{x[0]}" for x in result])
    return result


def get_random_list_elem(target_list: list):
    """Helper method that randomly chooses an element from a list. numpy.random.choice() doesn't 
    like to take elements from list()'s of tuples so this is the work around."""
    return target_list[np.random.randint(len(target_list))]  #np.random.choice only works on 1d


def pop_random_list_elem(target_list: list):
    """Helper method that randomly pop an element from a list. (delete from original list)
    numpy.random.choice() doesn't like to take elements from list()'s of tuples so
    this is the work around."""
    return target_list.pop(np.random.randint(len(target_list)))


def product_lists_allow_empty(list_of_lists: List[list]) -> List[list]:
    """product lists in a list and give all possible cases of one element from each.
    Different from itertool.product, not getting any element from a list is allowed.
    Note the result will include an empty list"""
    list_of_lists_copy = copy.deepcopy(list_of_lists)
    list_of_lists_copy = [each_list + [GHOST_LIST_ELEMENT] for each_list in list_of_lists_copy]
    return _product_lists_w_each_empty_ele(iter(list_of_lists_copy))


def _product_lists_w_each_empty_ele(list_of_lists: Iterable[list]) -> List[list]:
    """a sub-function used for product_list_allow_empty"""
    result_w_none = itertools.product(*list_of_lists)
    result = []
    count = 0
    for sublists in result_w_none:
        count += 1
        result_sublist = []
        for ele in sublists:
            if ele is GHOST_LIST_ELEMENT:
                continue
            result_sublist.append(ele)
        result.append(result_sublist)
    return result


def if_list_contain_repeating_element(target_list: list) -> bool:
    """check if the target list contains any repeating elements"""
    return len(target_list) != len(set(target_list))


def list_remove_adjacent_duplicates(target_list: list) -> list:
    """turn a new list removing all adjacent_duplicates. keep the 1st one.
    keep the order."""
    result = [target_list[0]]
    for i in range(1, len(target_list)):
        if target_list[i] != target_list[i - 1]:
            result.append(target_list[i])
    return result


def num_ele_2d(list_2d: List[list]) -> int:
    """counter the number of elements in a 2d list"""
    return sum(map(lambda x: len(x), list_2d))


# = Museum of Function =
# This is an old slow but insteresting function so it kept it here
# def _product_lists_w_each_empty_ele(list_of_lists: Iterable[list]) -> List[list]:
#     """a sub-function used for product_list_allow_empty"""
#     curr_list = next(list_of_lists, None)
#     if not curr_list:
#         return [[]]
#     next_list = _product_lists_w_each_empty_ele(list_of_lists)
#     return [[x] + y if x != GHOST_LIST_ELEMENT else y for x in curr_list
#             for y in next_list]


# == Dict related ==
def get_copy_of_deleted_dict(orig_dict: Dict, del_key) -> Dict:
    """
    get a copy of the orig_dict and delete an item in the copy base on the del_key
    the del_key can be a list of keys
    """
    dict_copy = copy.deepcopy(orig_dict)
    if isinstance(del_key, list):
        for single_key in del_key:
            del dict_copy[single_key]
    else:
        del dict_copy[del_key]

    return dict_copy


def swapped_dict(orig_dict: Dict) -> Dict:
    """get a swapped dictionary based on the original dictionary. The key and value are swapped
    in the new dictionary."""
    return {v : k for k, v in orig_dict.items()}


# == str related ==
def split_but_brackets(string: str, sep: str) -> List[str]:
    """split a string by {sep} but do not split anything in
    brackets. all {[( counts."""
    seperate_pattern = r"(?:[^"+ sep +"[{(]|(?:\{[^}]*\})|(?:\([^)]*\))|(?:\[[^\]]*\]))+[^"+ sep +r"]*"
    result = [i for i in re.findall(seperate_pattern, string)]
    return result


# == Class related ==
def get_str_for_print_class_var(cls) -> str:
    """return the str for printing out variables and values of cls to stdout"""
    result = ""
    class_members = dir(cls)
    # Filter out class variables (excluding methods and special members)
    class_variables = [member for member in class_members if not callable(getattr(cls, member)) and not member.startswith("__")]

    for var_name in class_variables:
        var_value = getattr(cls, var_name)
        result += f"{os.linesep}{var_name}: {os.linesep}{var_value}{os.linesep}"
    return result




# == context manager ==
class HiddenPrints:
    """block or redirect stdout/stderr prints to 'redirect'

    Example:
        with HiddenPrints(filename):
            print('this will go to filename')

    Note:
        This function will also redirect all the influenced logging handlers"""
    def __init__(self, redirect=os.devnull):
        self.redirect=redirect

    def __enter__(self):
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr
        self.redirect_stdout = open(self.redirect, 'w')
        self.redirect_stderr = open(self.redirect, 'w')
        sys.stdout = self.redirect_stdout
        sys.stderr = self.redirect_stderr
        # redirect any logging handler that lead to stdout/stderr
        all_loggers = [logging.getLogger("")] + [j for i, j in logging.Logger.manager.loggerDict.items()]
        for logger in all_loggers:
            if hasattr(logger, "handlers"):
                for handler in logger.handlers:
                    if hasattr(handler, "stream"):
                        if handler.stream is self.original_stdout:
                            handler.stream = self.redirect_stdout
                        if handler.stream is self.original_stderr:
                            handler.stream = self.redirect_stderr

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.redirect_stdout.close()
        self.redirect_stderr.close()
        sys.stdout = self.original_stdout
        sys.stderr = self.original_stderr
        # restore any logging handler that lead to redirect_stdout/stderr
        all_loggers = [logging.getLogger("")] + [j for i, j in logging.Logger.manager.loggerDict.items()]
        for logger in all_loggers:
            if hasattr(logger, "handlers"):
                for handler in logger.handlers:
                    if hasattr(handler, "stream"):
                        if handler.stream is self.redirect_stdout:
                            handler.stream = self.original_stdout
                        if handler.stream is self.redirect_stderr:
                            handler.stream = self.original_stderr


class EnablePropagate:
    """_LOGGER.propagate = True in the block"""
    def __init__(self, logger: logging.Logger) -> None:
        self.logger = logger

    def __enter__(self,):
        self.logger.propagate = True

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logger.propagate = False


class CaptureLogging:
    """redirect logging to a variable in the block"""
    def __init__(self, logger: logging.Logger) -> None:
        self.logger = logger

    def __enter__(self,):
        self.log_stream = StringIO()
        self.stream_handler = logging.StreamHandler(self.log_stream)
        self.old_handlers = self.logger.handlers[:]
        self.logger.handlers = [self.stream_handler]
        return self.log_stream

    def __exit__(self, exc_type, exc_val, exc_tb):
        # restore handlers
        self.logger.handlers = self.old_handlers
        # set position for read
        self.log_stream.seek(0)


# == misc ===
def timer(fn):
    """decodator for timing the run of the function {fn}"""

    def timer_inner(*args, **kwargs):
        start_time = time.perf_counter()
        to_execute = fn(*args, **kwargs)
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        _LOGGER.info("{0} took {1:.8f}s to execute".format(fn.__name__, execution_time))
        return to_execute

    return timer_inner


def get_localtime(time_stamp: float = None) -> str:
    """function that default return current locat time as formatted string.
    covert the {time_stamp} to localtime if provided"""
    if time_stamp is None:
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    else:
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time_stamp))


def get_itself(input_data: Any) -> Any:
    """a function that return the input itself"""
    return input_data


def save_obj(obj: Any, out_path: str):
    """save {obj} to the {out_path} as a .pickle file"""
    with open(out_path, "wb") as of:
        pickle.dump(obj, of)


def load_obj(in_path: str) -> Any:
    """load {obj} from the {in_path} as a .pickle file.
    Update: support loading a pickle file written multiple obj
    with 'ab'"""
    result = []
    with open(in_path, "rb") as f:
        while True:
            try:
                result.append(pickle.load(f))
            except EOFError:
                break

    if len(result) == 0:
        _LOGGER.error(f"The pickle file ({in_path}) have no content")
        raise EOFError
    elif len(result) == 1:
        return result[0]
    else:
        return result


def save_func_to_main(func: Callable, kwargs_file: str, out_path: str):
    """convert a function to a string of a python main script.
    This is commonly used when a function needs to be run by
    ARMer.
    Limitation: this func have to include import lines in itself.
    Will replace this with WorkUnit based solution later"""
    func_source = inspect.getsource(func)
    result = [
        "import sys",
        "from enzy_htp.core.general import load_obj",
        func_source,
        "if __name__ == '__main__':",
        f"    {func.__name__}(sys.argv, **load_obj('{kwargs_file}'))",
        "",
    ]
    write_lines(out_path, result)
