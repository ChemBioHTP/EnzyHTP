"""Testing the constants found in enzy_htp.chemical.metal

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-19
"""
import pytest
from typing import List

from enzy_htp.chemical import metal as mm

from util import all_caps



def test_variable_checks():
    """Basic testing for the metal mapper found in enzy_htp.chemical.metal"""
    mmapper_keys = list(mm.METAL_MAPPER.keys())
    mmapper_keys.remove("Na+")
    assert all_caps( mmapper_keys )
    assert max(list(map(len,mmapper_keys)))

    mcenter_keys = list(mm.METAL_CENTER_MAP.keys())
    assert all_caps( mcenter_keys )
    assert max(list(map(len,mcenter_keys)))
