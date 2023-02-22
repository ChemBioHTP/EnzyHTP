"""Testing the enzy_htp.molecular_mechanics.SystemConfig class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-03
"""
import enzy_htp
from enzy_htp._config import system_config
import pytest


def test_system_config_exists_in_singleton():
    """Checking that a system config is part of the singleton for the whole package."""
    assert enzy_htp.config._system


def test_assumed_attributes_and_values():
    """Checking that SystemConfig() has the assumed default attributes."""
    sconfig = enzy_htp._config.system_config.default_system_config()

    assert sconfig.N_CORES == 24
    assert sconfig.MEM_PER_CORE == 2000
    assert type(sconfig.WORK_DIR) == str and len(sconfig.WORK_DIR)
    assert sconfig.SCRATCH_DIR == (sconfig.WORK_DIR + "/scratch")


def test_getter_setter_round_trip():
    """Testing that the [] operator for getting and setting works via round-trip analysis."""
    sconfig = enzy_htp._config.system_config.default_system_config()

    sconfig['N_CORES'] = 1000
    assert sconfig['N_CORES'] == 1000

    sconfig['MEM_PER_CORE'] = 555
    assert sconfig['MEM_PER_CORE'] == 555

    sconfig['WORK_DIR'] = '/etc/'
    assert sconfig['WORK_DIR'] == '/etc/'

    sconfig['SCRATCH_DIR'] = '/sdir/'
    assert sconfig['SCRATCH_DIR'] == '/sdir/'


def test_errors_with_getters_and_setters():
    """Testing that the [] operator for getting and setting throws an error when appropriate."""

    sconfig = enzy_htp._config.system_config.default_system_config()
    with pytest.raises(SystemExit) as exe:
        sconfig['a.b'] = 5

    assert exe.type == SystemExit
    assert exe.value.code == 1

    with pytest.raises(SystemExit) as exe:
        temp_var = sconfig['a.b']

    assert exe.type == SystemExit
    assert exe.value.code == 1


def test_default_system_config_unique():
    """Checking that the config created by default_system_config() really is unique."""
    sysconfig1 = enzy_htp._config.system_config.default_system_config()
    assert id(sysconfig1) != id(enzy_htp.config._system)
    assert id(sysconfig1) != id(enzy_htp._config.system_config.default_system_config())
