"""Testing the enzy_htp.molecular_mechanics.AmberConfig class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-03
"""
import pytest

from pathlib import Path
import enzy_htp
from enzy_htp.core import file_system as fs
from enzy_htp._interface import AmberInterface
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


def test_invalid_box_type_throws():
    """Checking that the AmberConfig() setter throws when an invalid box type is set."""
    ac : mm.AmberConfig = mm.default_amber_config()
    with pytest.raises(TypeError) as exc:
        ac['BOX_TYPE'] = 'wrong'
    assert exc.type == TypeError

def test_valid_box_type():
    """Checking that the AmberConfig() setter doesn't throw for allowed box types."""
    ac : mm.AmberConfig = mm.default_amber_config()
    ac['BOX_TYPE'] = 'box'
    assert ac.valid_box_type() 
    ac['BOX_TYPE'] = 'oct'
    assert ac.valid_box_type() 

def test_required_executables():
    """Checking that the AmberConfig.required_executables() has the correct values."""
    ac : mm.AmberConfig = mm.default_amber_config()
    assert ac.required_executables() == ["$AMBERHOME/bin/sander.MPI","$AMBERHOME/bin/pmemd.cuda", "tleap", "ampdb", "parmchk2", "antechamber"]

def test_required_env_vars():
    """Checking that the AmberConfig.required_env_vars() contains the correct values."""
    ac : mm.AmberConfig = mm.default_amber_config()
    assert ac.required_env_vars() == ["AMBERHOME"]


def test_get_engine_valid():
    """Checking that AmberConfig.get_engine() works for valid inputs of "CPU" and "GPU"."""
    ac : mm.AmberConfig = mm.default_amber_config()
    assert ac.get_engine("CPU") == "$AMBERHOME/bin/sander.MPI"
    assert ac.get_engine("GPU") == "$AMBERHOME/bin/pmemd.cuda"

def test_get_engine_invalid():
    """Checking that AmberConfig.get_engine() fails for an invalid input"""
    ac : mm.AmberConfig = mm.default_amber_config()
    with pytest.raises(TypeError) as exc:
         _ = ac.get_engine('dne') 
    assert exc.type == TypeError

def test_round_trip_getting_and_setting():
    """Testing that bracket accession can be used in a round trip manner."""

    ac : mm.AmberConfig = mm.default_amber_config()
    ac['HOME'] = 'new_home'
    assert ac.HOME == 'new_home'
    assert ac['HOME'] == 'new_home'

