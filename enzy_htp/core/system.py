"""Module holds functions for manipulation of files, directories and other system related items. 

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-03-14
"""
import os, sys
import shutil
import datetime
from typing import List
from pathlib import Path
from .logger import _LOGGER


def safe_rm(fname: str) -> None:
    """Removes a file if and only if the directory already exists."""
    if os.path.exists(fname):
        os.remove(fname)


def safe_rmdir(dirname: str) -> None:
    """Removes a directory if and only if the directory already exists."""
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)


def safe_mkdir(dirname: str) -> None:
    """Makes a directory if and only if the directory does not already exist. Creates parents as needed"""
    if not os.path.isdir(dirname):
        Path(dirname).mkdir(parents=True, exist_ok=True)


def base_file_name(fname: str) -> str:
    """Given a filename, gets just the base name without extension or leading directories.
    Does not check if filex exists.
    """
    base_file = os.path.basename(fname)
    return os.path.splitext(base_file)[0]


def get_file_ext(fname: str) -> str:
    """TODO"""
    return Path(fname).suffix


def get_current_time() -> str:
    """Gets current system time in format YYYY_MM_DD_H_m"""
    return datetime.datetime.now().strftime('%Y_%m_%d_%H_%M')


def lines_from_file(fname: str) -> List[str]:
    """Extracts and returns lines from supplied filename. Returns empty list() if file odes not exist."""
    if not os.path.exists(fname):
        _LOGGER.error(f"The file {fname} does not exist.")
        return list()
    fh = open(fname, 'r')
    result = fh.read().splitlines()
    fh.close()
    return result


def write_lines(fname: str, lines: List[str]) -> None:
    """Writes lines to specified file, checking if file exists first and warning if it does."""
    if os.path.exists(fname):
        _LOGGER.warn(f"The file '{fname}' exists and will be overwritten")
    fh = open(fname, 'w')
    fh.write('\n'.join(lines))
    fh.close()


class HiddenPrints:
    """
    HiddenPrints(redirect=os.devnull)
    -----
    block or redirect stdout prints to "redirect"
    """

    def __init__(self, redirect=os.devnull):
        self.redirect = redirect

    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(self.redirect, "w")
        sys.stderr = open(self.redirect, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr
