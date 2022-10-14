"""Module holds functions for manipulation of files, directories and other system related items. 

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-03-14
"""
import os
import sys
import shutil
import datetime
from typing import List, Any, Dict
from pathlib import Path

from .logger import _LOGGER


def safe_rm(fname: str) -> None:
    """Removes a file if and only if the directory already exists."""
    # TODO(CJ): need to check if it is a directory. maybe re-direct if so?
    if os.path.exists(fname):
        os.remove(fname)


def is_empty_dir(dir_path: str) -> bool:
    """
    check if the dir_path is an empty dir
    """
    if os.path.isdir(dir_path):
        return not os.listdir(dir_path)
    else:
        _LOGGER.debug(f"No such directory: {dir_path}")
        return False


def safe_rmdir(dirname: str, empty_only: bool = False) -> None:
    """Removes a directory if and only if the directory already exists."""
    if os.path.isdir(dirname):
        if empty_only and not is_empty_dir(dirname):
            _LOGGER.debug(
                f"{dirname} is not empty. turn empty_only off to force remove it."
            )
        else:
            shutil.rmtree(dirname)


def safe_mkdir(dirname: str) -> None:
    """Makes a directory if and only if the directory does not already exist. Creates parents as needed."""
    if not os.path.isdir(dirname):
        Path(dirname).mkdir(parents=True, exist_ok=True)


def base_file_name(fname: str) -> str:
    """Given a filename, gets just the base name without extension or leading directories.
    Does not check if file exists.
    """
    base_file = os.path.basename(fname)
    return os.path.splitext(base_file)[0]


def remove_ext(fname: str) -> str:
    """Given a file name, returns the full filename EXCEPT the suffx."""
    return os.path.splitext(fname)[0]


def get_file_ext(fname: str) -> str:
    """Gets the file extension for a supplied file name."""
    return Path(fname).suffix


def get_current_time() -> str:
    """Gets current system time in format YYYY_MM_DD_H_m"""
    return datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")


def lines_from_file(fname: str) -> List[str]:
    """Extracts and returns lines from supplied filename. Returns empty list() if file odes not exist."""
    if not os.path.exists(fname):
        _LOGGER.error(f"The file {fname} does not exist.")
        return list()
    fh = open(fname, "r")
    result = fh.read().splitlines()
    fh.close()
    return result


def write_lines(fname: str, lines: List[str]) -> None:
    """Writes lines to specified file, checking if file exists first and warning if it does. Assumes no newlines."""
    # TODO(CJ) check if binary file and dont return if so
    if os.path.exists(fname):
        _LOGGER.warning(f"The file '{fname}' exists and will be overwritten")
    fh = open(fname, "w")
    fh.write("\n".join(lines))
    fh.close()


def write_data(outfile: str, tag: Any, data: Dict) -> str:
    #TODO(CJ): add the doc-string and also unittests
    tag: str = repr(tag)
    fh = open(outfile, 'a')
    fh.write("===TAG===\n")
    fh.write(f"{tag}\n")
    for key, value in data.items():
        fh.write(f"---{key}---\n")
        fh.write(f"{repr(value)}\n")
    fh.close()

    return outfile


def get_valid_temp_name(fname: str) -> None:
    """find a vaild name for a temp file"""
    idx = 0
    suffix = ''.join(Path(fname).suffixes)
    result_fname = fname
    while os.path.isfile(result_fname):
        idx += 1
        result_fname = f"{fname[:-len(suffix)]}_{str(idx)}{suffix}"
    return result_fname
