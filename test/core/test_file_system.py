"""Testing the utility functions in the enzy_htp.core.system sub module.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os

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


# TODO(CJ) add tests for remove_ext
