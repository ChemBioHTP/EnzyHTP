"""Testing the utility functions in the enzy_htp.core.system sub module.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import multiprocessing
import os
import time
from pathlib import Path

import enzy_htp.core.file_system as fs


def make_test_file(fname: str, contents: str) -> None:
    """Helper function that writes supplied contents to a supplied filename."""
    fh = open(fname, "w")
    fh.write(contents)
    fh.close()


CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)


def test_safe_rm():
    """Ensuring the safe_rm() method removes files."""
    test_fname = f"{CURR_DIR}/__test.txt"
    make_test_file(test_fname, "test")

    assert os.path.exists(test_fname)
    assert not fs.safe_rm(test_fname)
    assert not os.path.exists(test_fname)

    assert not fs.safe_rm(test_fname)


def test_safe_rmdir():
    """Making sure that the safe_rmdir() method works for a variety of situations."""

    test_dir = f"{CURR_DIR}/__test_dir/"
    test_file = f"{test_dir}/test.txt"
    test_dir_child = f"{test_dir}/child1/"
    test_file_child = f"{test_dir_child}/test.txt"

    os.mkdir(test_dir)
    os.mkdir(test_dir_child)
    make_test_file(test_file, "test")
    make_test_file(test_file_child, "test")

    assert os.path.isdir(test_dir)
    assert os.path.isdir(test_dir_child)
    assert os.path.exists(test_file)
    assert os.path.exists(test_file_child)

    fs.safe_rmdir(test_dir_child)

    assert os.path.isdir(test_dir)
    assert not os.path.isdir(test_dir_child)
    assert os.path.exists(test_file)
    assert not os.path.exists(test_file_child)

    os.mkdir(test_dir_child)
    make_test_file(test_file_child, "test")

    fs.safe_rmdir(test_dir)

    assert not os.path.isdir(test_dir)
    assert not os.path.isdir(test_dir_child)
    assert not os.path.exists(test_file)
    assert not os.path.exists(test_file_child)


def test_safe_mkdir():
    """Ensuring that the safe_mkdir() function works for simple and complex cases."""

    test_dir = f"{CURR_DIR}/__test_dir/"
    test_dir_child1 = f"{CURR_DIR}/__test_dir/_child1/"
    test_dir_child2 = f"{CURR_DIR}/__test_dir/_child1/_child1"

    assert not os.path.isdir(test_dir)
    assert not os.path.isdir(test_dir_child1)
    assert not os.path.isdir(test_dir_child2)

    fs.safe_mkdir(test_dir_child2)

    assert not fs.safe_mkdir(test_dir_child2)

    assert os.path.isdir(test_dir)
    assert os.path.isdir(test_dir_child1)
    assert os.path.isdir(test_dir_child2)

    fs.safe_rmdir(test_dir)


def test_base_file_name():
    """Checking that base_file_name() works for a few edge cases."""

    assert fs.base_file_name("file1.txt") == "file1"
    assert fs.base_file_name("/path/to/file/file1.txt") == "file1"

    assert fs.base_file_name(".file.txt") != fs.base_file_name("file.txt")


def test_get_file_ext():
    """Ensure that get_file_ext() works for a few edge cases."""
    assert fs.get_file_ext("file1.txt") == ".txt"
    assert fs.get_file_ext(".bashrc") == ""
    assert fs.get_file_ext("/path/to/file1.txt") == ".txt"
    assert fs.get_file_ext("/path/to/file1.txt.bk") == ".bk"
    assert fs.get_file_ext("pdb.PDB") != fs.get_file_ext("pdb.pdb")


def test_file_name_round_trip():
    """Checking that the file name round trip works for base_file_name() and get_file_ext()"""
    fname1 = "test.txt"
    fname2 = "/path/to/test.txt"
    assert fs.base_file_name(fname1) + fs.get_file_ext(fname1) == fname1
    assert fs.base_file_name(fname2) + fs.get_file_ext(fname2) == fname1


def test_get_current_time():
    """Hard to test, but just checking that get_curren_time() works."""
    curr_time = fs.get_current_time()
    assert curr_time.count("_") == 4
    assert len(curr_time) == 16
    lens = list(map(len, curr_time.split("_")))
    assert lens == [4, 2, 2, 2, 2]


def test_lines_from_file():
    """Checks that lines_from_file() works properly."""
    dne = "__doesnotest.txt"
    assert not len(fs.lines_from_file(dne))
    lines = ["line1", "line2", "line3"]
    fname1 = f"{CURR_DIR}/test_file.txt"
    make_test_file(fname1, "\n".join(lines))
    assert fs.lines_from_file(fname1) == lines
    fs.safe_rm(fname1)
    assert not os.path.exists(fname1)


def test_write_lines():
    """Checking that test_write_lines() functions correctly."""
    contents = ["line1", "line2", "line3"]
    fname1 = f"{CURR_DIR}/test_fle.txt"
    make_test_file(fname1, "")
    assert fs.lines_from_file(fname1) != contents
    fs.write_lines(fname1, contents)
    assert fs.lines_from_file(fname1) == contents
    fs.safe_rm(fname1)
    assert not os.path.exists(fname1)


def test_get_valid_temp_name():
    fname1 = f"{CURR_DIR}/test.test.test.py"
    make_test_file(fname1, "")
    assert fs.get_valid_temp_name(fname1) == f"{CURR_DIR}/test_000001.test.test.py"
    fs.safe_rm(fname1)
    assert not os.path.exists(fname1)

def test_get_valid_temp_name_symlink():
    sym_fname1 = f"{CURR_DIR}/test.test.test.py.lnk"
    fname1 = f"{CURR_DIR}/test.test.test.py"
    Path(sym_fname1).symlink_to(fname1)
    
    assert fs.get_valid_temp_name(sym_fname1, is_symlink=True) == f"{CURR_DIR}/test_000001.test.test.py.lnk"
    fs.safe_rm(sym_fname1)
    assert not os.path.lexists(sym_fname1)

def test_get_valid_temp_name_dir():
    fname1 = f"{CURR_DIR}/test_get_valid_temp_name/"
    fname2 = f"{CURR_DIR}/test_get_valid_temp_name_000001/////"
    fs.safe_mkdir(fname1)
    assert fs.get_valid_temp_name(fname1) == f"{CURR_DIR}/test_get_valid_temp_name_000001"
    fs.safe_mkdir(fname2)
    assert fs.get_valid_temp_name(fname1) == f"{CURR_DIR}/test_get_valid_temp_name_000002"
    fs.safe_rmdir(fname1)
    fs.safe_rmdir(fname2)
    assert not os.path.exists(fname1)
    assert not os.path.exists(fname2)

def test_is_path_exist_dir():
    fname1 = f"{CURR_DIR}/test_get_valid_temp_name_real/"
    fname2 = f"{CURR_DIR}/test_get_valid_temp_name_fake/"
    fs.safe_mkdir(fname1)
    assert fs.is_path_exist(fname1)
    assert not fs.is_path_exist(fname2)
    fs.safe_rmdir(fname1)
    assert not os.path.exists(fname1)
    assert not os.path.exists(fname2)

def test_is_path_exist_symlink():
    sym_fname1 = f"{CURR_DIR}/test.test.fake.py.lnk"
    fname1 = f"{CURR_DIR}/test.test.fake.py"

    Path(sym_fname1).symlink_to(fname1)
    assert fs.is_path_exist(sym_fname1, is_symlink=True)
    assert not fs.is_path_exist(sym_fname1, is_symlink=False)
    fs.safe_rm(sym_fname1)
    assert not os.path.lexists(sym_fname1)

def test_is_path_exist():
    fname1 = f"{CURR_DIR}/test.test.test.gjf"
    fname2 = f"{CURR_DIR}/test.test.fake.gjf"
    
    make_test_file(fname1, "")
    assert fs.is_path_exist(fname1)
    assert not fs.is_path_exist(fname2)
    fs.safe_rm(fname1)
    assert not os.path.exists(fname1)

def test_is_path_exist_multi_ext():
    fname1 = f"{CURR_DIR}/test.gjf"
    fname2 = f"{CURR_DIR}/test.out"
    
    make_test_file(fname2, "")
    assert fs.is_path_exist(fname1, ext_set=[".chk", ".out"])
    fs.safe_rm(fname2)
    assert not os.path.exists(fname2)

def test_clean_temp_file_n_dir():
    temp_dir_path = f"{CURR_DIR}/temp/"
    temp_path_list = [temp_dir_path]
    fs.safe_mkdir(temp_dir_path)
    for i in range(3):
        temp_file = fs.get_valid_temp_name(f"{temp_dir_path}/temp.txt")
        with open(temp_file, "w") as of:
            of.write("test")
        temp_path_list.append(temp_file)
    for temp_path in temp_path_list:
        assert os.path.exists(temp_path)
    fs.clean_temp_file_n_dir(temp_path_list)
    for temp_path in temp_path_list:
        assert not os.path.exists(temp_path)


def test_all_file_in_dir():
    """test the function using dummy folder and files"""
    test_dirname = f"{CURR_DIR}/data/dummy_dir"
    recursive_files = fs.all_file_in_dir(test_dirname, recursive=True)
    assert len(recursive_files) == 2

    nonrecursive_files = fs.all_file_in_dir(test_dirname, recursive=False)
    assert len(nonrecursive_files) == 1


def test_lock():
    """test the homemade file lock"""
    def worker_1():
        f = open("test.file", "wb")
        fs.lock(f)
        f.write(b"Hi")
        time.sleep(5)
        fs.unlock(f)
        f.close()

    def worker_2(result_queue):
        f = open("test.file", "rb")
        while fs.is_locked(f):
            print(1)
            time.sleep(0.5)
        content = f.read()
        result_queue.put(content)
        f.close()

    result_queue = multiprocessing.Queue()
    process1 = multiprocessing.Process(target=worker_1)
    process2 = multiprocessing.Process(target=worker_2, args=(result_queue,))
    process1.start()
    time.sleep(0.1)
    process2.start()
    process1.join()
    process2.join()
    fs.safe_rm("test.file")

    assert result_queue.get() == b"Hi"

# TODO(CJ) add tests for remove_ext
