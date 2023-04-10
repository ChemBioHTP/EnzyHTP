"""Testing the enzy_htp._config.AmberConfig class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-03
"""

from pathlib import Path
from enzy_htp._config import amber_config
import pytest
import enzy_htp._config.amber_config as ac


def test_default_amber_config():
    """Ensuring the default_amber_config() method returns deepcopied AmberConfig() objects."""
    ac1: ac.AmberConfig = ac.default_amber_config()
    ac2: ac.AmberConfig = ac.default_amber_config()
    assert ac1
    assert ac2
    assert id(ac1) != id(ac2)


def test_required_attributes():
    """Checking that required attributes for the AmberConfig() object exist."""
    ac1: ac.AmberConfig = ac.default_amber_config()
    assert ac1.HOME
    assert ac1.CPU_ENGINE
    assert ac1.GPU_ENGINE
    assert ac1.BOX_TYPE
    assert ac1.BOX_SIZE
    assert ac1.CONF_MIN
    assert ac1.CONF_HEAT
    assert ac1.CONF_EQUI
    assert ac1.CONF_PROD


def test_invalid_box_type_throws():
    """Checking that the AmberConfig() setter throws when an invalid box type is set."""
    ac1: ac.AmberConfig = ac.default_amber_config()
    with pytest.raises(TypeError) as exc:
        ac1['BOX_TYPE'] = 'wrong'
    assert exc.type == TypeError


def test_valid_box_type():
    """Checking that the AmberConfig() setter doesn't throw for allowed box types."""
    ac1: ac.AmberConfig = ac.default_amber_config()
    ac1['BOX_TYPE'] = 'box'
    assert ac1.valid_box_type()
    ac1['BOX_TYPE'] = 'oct'
    assert ac1.valid_box_type()


def test_required_executables():
    """Checking that the AmberConfig.required_executables() has the correct values."""
    ac1: ac.AmberConfig = ac.default_amber_config()
    assert ac1.required_executables() == [
        "sander", "pmemd.cuda", "tleap", "ambpdb", "parmchk2", "antechamber", "cpptraj"
    ]


def test_required_env_vars():
    """Checking that the AmberConfig.required_env_vars() contains the correct values."""
    ac1: ac.AmberConfig = ac.default_amber_config()
    assert ac1.required_env_vars() == ["AMBERHOME"]


def test_get_engine_valid():
    """Checking that AmberConfig.get_engine() works for valid inputs of "CPU" and "GPU"."""
    ac1: ac.AmberConfig = ac.default_amber_config()
    assert ac1.get_engine("CPU") == "sander"
    assert ac1.get_engine("GPU") == "pmemd.cuda"


def test_get_engine_invalid():
    """Checking that AmberConfig.get_engine() fails for an invalid input"""
    ac1: ac.AmberConfig = ac.default_amber_config()
    with pytest.raises(TypeError) as exc:
        _ = ac1.get_engine('dne')
    assert exc.type == TypeError


def test_round_trip_getting_and_setting():
    """Testing that bracket accession can be used in a round trip manner."""

    ac1: ac.AmberConfig = ac.default_amber_config()
    ac1['HOME'] = 'new_home'
    assert ac1.HOME == 'new_home'
    assert ac1['HOME'] == 'new_home'
