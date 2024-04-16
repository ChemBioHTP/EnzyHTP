"""Testing the enzy_htp._config.RosettaConfig class.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-03-04
"""

import pytest
import os

from enzy_htp import config
rc = config.rosetta


def test_get_cart_ddg_exe():
    """Ensuring the get_cart_ddg_exe() method returns the correct path."""
    os.environ["ROSETTA3"] = "/data/yang_lab/Common_Software/Rosetta3.9/main"
    result = rc.get_cart_ddg_exe()
    assert result == "/data/yang_lab/Common_Software/Rosetta3.9/main/source/bin/cartesian_ddg.static.linuxgccrelease"
