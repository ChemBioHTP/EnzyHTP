"""Testing the enzy_htp.molecular_mechanics.AmberConfig class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-03
"""
from pathlib import Path
from enzy_htp.core import file_system as fs
from enzy_htp import AmberInterface
from enzy_htp import molecular_mechanics as mm
from pathlib import Path


def test_default_amber_config():
    """Ensuring the default_amber_config() method returns deepcopied AmberConfig() objects."""
    ac1 : mm.AmberConfig = mm.default_amber_config()
    ac2 : mm.AmberConfig = mm.default_amber_config()
    assert ac1
    assert ac2
    assert id(ac1) != id(ac2)

def test_required_attributes():
    """Checking that required attributes for the AmberConfig() object exist."""
    ac1 : mm.AmberConfig = mm.default_amber_config()
    assert ac1.HOME
    assert ac1.CPU_ENGINE
    assert ac1.GPU_ENGINE
    assert ac1.BOX_TYPE
    assert ac1.BOX_SIZE
    assert ac1.CONF_MIN
    assert ac1.CONF_HEAT
    assert ac1.CONF_EQUI
    assert ac1.CONF_PROD

