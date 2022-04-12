"""Testing the functions found in the enzy_htp.core.settings submodule.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-11
"""
import os,sys
import pytest
from enzy_htp.core import settings as ss

def test_get_version():
    """Checking that the version number is correct."""
    assert ss.get_version() == ss.EnzyHTPVersion(major=0,minor=1,patch=0)


def test_get_version_str():
    """Checking that the version number is correct in the string format.."""
    assert ss.version_to_str(ss.get_version()) == "0.1.0" 

def test_version_to_str_truncation():
	"""Checking that the version_to_str() method truncates correctly."""

	assert ss.version_to_str(ss.EnzyHTPVersion(major=1.5,minor=0,patch=0)) == "1.0.0"
	assert ss.version_to_str(ss.EnzyHTPVersion(major=1.5,minor=2.9,patch=0)) == "1.2.0"
	assert ss.version_to_str(ss.EnzyHTPVersion(major=1.5,minor=2.9,patch=3.2)) == "1.2.3"


def test_is_compatible_os():
    """Checking that the is_compatible_os() method works."""
    cached_platform : str = sys.platform
    sys.platform = "linux"
    assert ss.is_compatible_os() == True
    sys.platform = "darwin"
    assert ss.is_compatible_os() == True
    sys.platform = "windows"
    assert ss.is_compatible_os() == False
    sys.platform = "win32"
    assert ss.is_compatible_os() == False
    sys.platform = "NA"
    assert ss.is_compatible_os() == False
    sys.platform = cached_platform


def test_data_dir():
    """Ensuring that the data_dir() method works."""
    result = ss.data_dir()
    tks = list(filter(len,result.split('/')))
    assert tks[-1] == 'data'
    assert tks[-1] == 'enzy_htp'
