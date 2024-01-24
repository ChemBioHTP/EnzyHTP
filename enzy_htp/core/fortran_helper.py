"""this module contain helper functions for stuffs related to Fortran language.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Author: QZ Shao, <shaoqz@icloud.com>
Date: 2024-01-23
"""
from typing import Any, Tuple, List
import re

from .logger import _LOGGER

def parse_data(fmt: str, data: str, reduce: bool = False) -> List[List]:
    """helper function that parses a data section from a Fortran-based data file.
    (e.g.: the prmtop file from Amber MD)

    Args:
        fmt:
            the FORMAT() string. (e.g.: "FORMAT(I5, 2F10.2, A10, /, T10, L5)")
        data:
            the data content.
        reduce:
            whether reduce the layer of the resulting data list.
            i.e.: if there is only 1 element in the data point list, the list will be reduced
            to that element.

    Returns:
        a list of tuple that each contain a data point following the format.
    """
    result = []
    
    format_list = parse_format(fmt)
    total_width = get_total_width(format_list)

    raw_data_list: List[str] = list()
    for idx in range(0, len(data), total_width):
        raw_data_list.append(data[idx:idx + total_width])

    for raw_data in raw_data_list:
        idx = 0
        result_data = []
        for repeat, type_ctor, width, digit in format_list:
            while repeat > 0:
                if idx + width > len(raw_data):
                    _LOGGER.info("repeat not fullfilled in data. make sure it is what you want")
                    break
                raw_data_ele = raw_data[idx: idx + width].strip()
                result_data.append(type_ctor(raw_data_ele))
                repeat -= 1
                idx += width
        if reduce and len(result_data) == 1:
            result_data = result_data[0]
        result.append(result_data)
    
    return result

def parse_format(fmt: str) -> List[Tuple[Any, int]]:
    """helper function that parses a FORMAT() notation from a Fortran-based data file.
    (e.g.: the prmtop file from Amber MD)
    Identifies (repeat, type, width, digit) of the format. Supplied string can contain whitespace.
    In the case that the supplied format string does not contain the pattern 'FORMAT(<>)',
    the result (None, -1) is returned.

    Args:
        fmt: the FORMAT() string. (e.g.: "FORMAT(I5, 2F10.2, A10, /, T10, L5)")

    Returns:
        a list of tuple that each contain
        (number of repeat, type ctor, data width, data digit).
    """
    # san check
    if (not fmt.strip().startswith('FORMAT(')) or (not fmt.strip().endswith(')')):
        return (-1, None, -1, -1)

    # clean up
    fmt = fmt.strip().replace('FORMAT(', '').replace(')', '')

    result = []

    fmt_list = [i.strip() for i in fmt.split(",")]
    for fmt in fmt_list:
        result.append(parse_single_format(fmt))
    
    return result

def parse_single_format(fmt: str) -> Tuple[Any, int]:
    """parse a single format notation such as "I5"
    return a tuple() of length 4 containing the
    (number of repeat, type ctor, data width, data digit)."""
    # limitation check
    unsupported_types = [
        "T", "/", ":"
    ]
    for type_chr in unsupported_types:
        if type_chr in fmt:
            _LOGGER.error(f"dont support parsing '{type_chr}' yet (input: {fmt})")
            raise ValueError
    
    fortran_type_mapper = {
        "I": int,       # Integer
        "F": float,     # Floating-point number
        "D": float,     # Double-precision floating-point number (in Python, float is typically double-precision)
        "C": complex,   # Complex number
        "L": bool,      # Logical type (Boolean)
        "A": str,       # Character string
        "E": float,     # Exponential format floating-point number
        "G": float      # General floating-point number (use float in Python)
    }
    """This dictionary maps Fortran type characters to their corresponding Python type constructors."""

    pattern = r"(\d*)(\w)(\d+)\.?(\d*)"
    match_result = re.match(pattern, fmt)
    if match_result:
        repeat, type_chr, width, digit = match_result.groups()
    else:
        _LOGGER.error(f"{fmt} dont follow pattern {pattern}. parse failed.")
        raise ValueError

    type_chr = type_chr.upper()
    type_ctor = fortran_type_mapper[type_chr]
    if not repeat:
        repeat = 1
    if not digit:
        digit = 0

    return int(repeat), type_ctor, int(width), int(digit)

def get_total_width(fmt_list: List[Tuple]) -> int:
    """determine the total width of a datapoint composed by
    a list of data elements that each have a width.
    the total width is the sum of the width of each element"""
    result = 0
    for repeat, ctor, width, digit in fmt_list:
        while repeat > 0:
            result += width
            repeat -= 1
    return result
