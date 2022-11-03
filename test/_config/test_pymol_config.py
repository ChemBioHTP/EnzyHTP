"""Testing the enzy_htp._config.PyMOLConfig class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-11-03
"""

import pytest
from enzy_htp._config import pymol_config 
import enzy_htp._config.pymol_config as pc


def test_default_pymol_config():
    """Ensuring the default_pymol_config() method returns deepcopies PyMOLConfig() objects."""   
    pc1 : pc.PyMOLConfig = pc.default_pymol_config()
    pc2 : pc.PyMOLConfig = pc.default_pymol_config()
    assert pc1
    assert pc2
    assert id(pc1) != id(pc2)

