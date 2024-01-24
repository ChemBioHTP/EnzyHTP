"""Testing enzy_hpt.core.fortran_helper.py

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-01-23
"""
import enzy_htp.core.fortran_helper as fh

def test_parse_single_format():
    """as name"""
    test = ""
    assert fh.parse_single_format(test) == ()

# TODO
def test_parse_format_uppercase():
    """Testing that the parse_format() method works correctly for uppercase inputs."""
    assert fh.parse_format('%FORMAT(20A4)') == (str, 4)

    assert fh.parse_format('%FORMAT(5E16.8)') == (float, 16)

    assert fh.parse_format('%FORMAT(10I8)') == (int, 8)


def test_parse_format_lowercase():
    """Testing that the parse_format() method works correctly for lowercase inputs."""
    assert fh.parse_format('%FORMAT(20a4)') == (str, 4)

    assert fh.parse_format('%FORMAT(5e16.8)') == (float, 16)

    assert fh.parse_format('%FORMAT(10i8)') == (int, 8)


def test_parse_format_bad_input():
    """Testing that the parse_format() method returns the correct (None,-1) for bad inputs."""


    assert fh.parse_format('') == (None, -1)

    assert fh.parse_format('20a4') == (None, -1)

    assert fh.parse_format('FORMAT(20a4)') == (None, -1)

    assert fh.parse_format('%FORMAT(20a4') == (None, -1)

