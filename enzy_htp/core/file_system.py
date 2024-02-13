"""Module holds functions for manipulation of files, directories and other system related items. 

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Author: QZ Shao, <shaoqz@icloud.com>
Date: 2022-03-14
"""
from io import IOBase
import os
import glob
import shutil
import datetime
from typing import List, Any, Dict
from pathlib import Path

from .logger import _LOGGER

# == general ==
def safe_rm(fname: str) -> None:
    """Removes a file if and only if the directory already exists. Provides a warning if the 
    supplied path is a directory."""
    if os.path.isdir(fname):
        _LOGGER.warning(f"The supplied path '{fname}' is a directory and cannot be removed by enzy_htp.core.file_system.safe_rm()")
        return

    if os.path.lexists(fname):
        os.remove(fname)

def safe_mv(src: str, dest: str) -> str:
    """Robust move method that replicates shell 'mv' command. Creates temporary directories
    when needed.
    
    Args:
        src: Initial path to file object.
        dest: Where the src will be placed. Either a new filename or destination directory.

    Returns:
        Path to the final location of the src file object.        
    """
    src = Path(src)
    dest = Path(dest)

    if dest.is_dir():
        safe_mkdir(dest)
        dest = dest / src.name

    return str(shutil.move(src, dest))

# == dir ==
def is_empty_dir(dir_path: str) -> bool:
    """Is the supplied directory empty?"""
    if os.path.isdir(dir_path):
        return not os.listdir(dir_path)
    else:
        _LOGGER.debug(f"No such directory: {dir_path}")
        return False


def safe_rmdir(dirname: str, empty_only: bool = False) -> None:
    """Removes a directory if and only if the directory already exists."""
    if os.path.isdir(dirname):
        if empty_only and not is_empty_dir(dirname):
            _LOGGER.debug(f"{dirname} is not empty. turn empty_only off to force remove it.")
        else:
            shutil.rmtree(dirname)


def safe_mkdir(dirname: str) -> None:
    """Makes a directory if and only if the directory does not already exist. Creates parents as needed."""
    if not os.path.isdir(dirname):
        Path(dirname).mkdir(parents=True, exist_ok=True)


def all_file_in_dir(dirname: str, recursive: bool=True):
    """get path of all files under a dir. filter out all subdir paths themselves"""
    if recursive:
        file_paths = glob.glob(f"{dirname}/**/*", recursive=True)
    else:
        file_paths = glob.glob(f"{dirname}/*")
    # filter out all dir
    file_paths = [f for f in file_paths if not os.path.isdir(f)]

    return file_paths

# == path/file ==
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


def lines_from_file(fname: str) -> List[str]:
    """Extracts and returns lines from supplied filename. Returns empty list() if file does not exist."""
    if not os.path.exists(fname):
        _LOGGER.error(f"The file {fname} does not exist.")
        return list()
    fh = open(fname, "r")
    result = fh.read().splitlines()
    fh.close()
    return result


def content_from_file(fname: str) -> str:
    """Extracts and returns the content from supplied filename. Returns empty str() if file does not exist."""
    #TODO(CJ): make unit tests for this
    if not os.path.exists(fname):
        _LOGGER.error(f"The file {fname} does not exist.")
        return str()
    fh = open(fname, "r")
    result = fh.read()
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


def get_valid_temp_name(fname: str, is_symlink: bool = False) -> str:
    """find a vaild name for a temporary file of {fname}. 
    If {fname} exists, add a index after it.
    Args:
        fname: the filename of need.
        is_symlink: whether the target filename is a symlink
    Return:
        result_fname: the valid filename that is unique to use"""
    idx = 0
    suffix = ''.join(Path(fname).suffixes)
    result_fname = fname
    if is_symlink:
        while Path(result_fname).is_symlink():
            idx += 1
            result_fname = f"{fname[:-len(suffix)]}_{str(idx)}{suffix}"
    else:
        while os.path.isfile(result_fname):
            idx += 1
            result_fname = f"{fname[:-len(suffix)]}_{str(idx)}{suffix}"
    return result_fname


def check_file_exists(fname: str, exit_script: bool = True) -> None:
    """Function that checks if a file exists. Will either exit the script or raise
    an error depending on the specified behavior.

    Args:
        fname: The str() name of the file to check.
        exit_script: A bool() specifying if the script should exit or raise a FileNotFoundError().
        
    Returns:
        Nothing.

    Raises:
        FileNotFoundError if exit_script is set to False and the file does not exist.
    """

    if Path(fname).exists():
        return

    if exit_script:
        _LOGGER.error(f"The file '{fname}' does not exist. Exiting...")
        exit(1)
    else:
        _LOGGER.error(f"The file '{fname}' does not exist. Exiting...")
        raise FileNotFoundError(f"The file '{fname}' does not exist.")


def check_not_empty(fname: str) -> None:
    """Function that checks if a file exists and is not empty. Will exit the script if either
    of the two conditions are false.
    
    Args:
        fname: The str() name of the file to check.

    Returns:
        Nothing.
    """

    check_file_exists(fname)

    temp = Path(fname)
    if temp.stat().st_size > 0:
        return

    _LOGGER.error(f"The file '{fname}' is empty. Exiting...")

    exit(1)


def has_content(fname: str) -> bool:
    """Method that checks if the supplied path exists and contains content. Returns bool with result

    Args:
        fname: The name of the file to check as a str().

    Returns:
        Whether the file has content.
    """
    fpath = Path(fname)

    return fpath.exists() and fpath.stat().st_size > 0


def clean_temp_file_n_dir(temp_path_list: List[str]) -> None:
    """clean temporary files created by EnzyHTP functions.
    removes file first and then dirs.
    Args:
        temp_path_list:
            a list of temp_path that could be either file or dir"""
    if _LOGGER.level > 10:  # not DEBUG or below
        dir_list = filter(os.path.isdir, temp_path_list)
        file_list = filter(os.path.isfile, temp_path_list)
        for file_path in file_list:
            safe_rm(file_path)
        for dir_path in dir_list:
            safe_rmdir(dir_path, empty_only=True)


# make own lock function that python ones not really work
def is_locked(f: IOBase) -> bool:
    """check whether the file is locked by another process"""
    lock_path = f"{os.path.abspath(f.name)}.lock"
    return Path(lock_path).exists()

def lock(f: IOBase):
    """lock a file by making a lock file
    * this does not work when the read happens too fast"""
    lock_path = f"{os.path.abspath(f.name)}.lock"
    if is_locked(f):
        _LOGGER.error(f"lock already exist in {lock_path}")
        raise IOError
    with open(lock_path, "w"):
        pass

def unlock(f: IOBase):
    """unlock a file by deleting a lock file"""
    lock_path = f"{os.path.abspath(f.name)}.lock"
    if is_locked(f):
        safe_rm(lock_path)
    else:
        _LOGGER.warning(f"unlocking non-existing lock: {lock_path}")

# == time ==
def get_current_time() -> str:
    """Gets current system time in format YYYY_MM_DD_H_m"""
    return datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")

