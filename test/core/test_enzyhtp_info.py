"""Testing enzy_hpt.core.enzyhtp_info.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-20
"""

from enzy_htp.core import enzyhtp_info

def test_enzyhtpversion_str():
    """test for the __str__ method of EnzyHTPVersion"""
    test_enzyhtp_version = enzyhtp_info.EnzyHTPVersion(0, 0, 0)
    assert str(test_enzyhtp_version) == "0.0.0"

def test_get_version():
    """test for get_version"""
    current_version = enzyhtp_info.VERSION
    # test by creating a pesudo testing version
    test_enzyhtp_version = enzyhtp_info.EnzyHTPVersion(0, 0, 0)
    enzyhtp_info.VERSION = test_enzyhtp_version
    assert enzyhtp_info.get_version() is test_enzyhtp_version
    # restore current version
    enzyhtp_info.VERSION = current_version
