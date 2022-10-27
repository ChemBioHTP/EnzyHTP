"""Testing the logger from enzy_hpt.core.logger.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os
import pytest
import logging

from enzy_htp.core import file_system as fs
from enzy_htp.core import logger as lg

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_FILE = f"{CURR_DIR}/__log.test"
fs.safe_rm(LOG_FILE)


def test_get_logger():
    """Testing that the logger can be acquired."""
    _LOG_ALIAS = None
    assert not os.path.exists(LOG_FILE)
    _LOG_ALIAS = lg.init_logger("unit-testing", LOG_FILE, start=True)
    assert os.path.exists(LOG_FILE)
    assert _LOG_ALIAS
    fs.safe_rm(LOG_FILE)
    assert not os.path.exists(LOG_FILE)


def test_log_data_debug_level():
    """Making sure that the logger works at debug level."""
    _LOG_ALIAS = lg.init_logger("unit-testing", LOG_FILE, start=True)
    assert not fs.lines_from_file(LOG_FILE)
    _LOG_ALIAS.setLevel(logging.DEBUG)
    _LOG_ALIAS.debug("testing the logger, debug")
    _LOG_ALIAS.info("testing the logger, info")
    _LOG_ALIAS.warning("testing the logger, warning")
    _LOG_ALIAS.error("testing the logger, error")
    _LOG_ALIAS.critical("testing the logger, info")

    lines = fs.lines_from_file(LOG_FILE)

    assert len(lines) == 5
    assert lines[0].find("DEBUG") != -1
    assert lines[1].find("INFO") != -1
    assert lines[2].find("WARNING") != -1
    assert lines[3].find("ERROR") != -1
    assert lines[4].find("CRITICAL") != -1

    fs.safe_rm(LOG_FILE)
    assert not os.path.exists(LOG_FILE)


def test_log_data_info_level():
    """Making sure that the logger works at info level."""
    _LOG_ALIAS = lg.init_logger("unit-testing", LOG_FILE, start=True)
    assert not fs.lines_from_file(LOG_FILE)
    _LOG_ALIAS.setLevel(logging.INFO)
    _LOG_ALIAS.debug("testing the logger, debug")
    _LOG_ALIAS.info("testing the logger, info")
    _LOG_ALIAS.warning("testing the logger, warning")
    _LOG_ALIAS.error("testing the logger, error")
    _LOG_ALIAS.critical("testing the logger, info")

    lines = fs.lines_from_file(LOG_FILE)

    assert len(lines) == 4
    assert lines[0].find("INFO") != -1
    assert lines[1].find("WARNING") != -1
    assert lines[2].find("ERROR") != -1
    assert lines[3].find("CRITICAL") != -1

    fs.safe_rm(LOG_FILE)
    assert not os.path.exists(LOG_FILE)


def test_log_data_warning_level():
    """Making sure that the logger works at warning level."""
    _LOG_ALIAS = lg.init_logger("unit-testing", LOG_FILE, start=True)
    assert not fs.lines_from_file(LOG_FILE)
    _LOG_ALIAS.setLevel(logging.WARNING)
    _LOG_ALIAS.debug("testing the logger, debug")
    _LOG_ALIAS.info("testing the logger, info")
    _LOG_ALIAS.warning("testing the logger, warning")
    _LOG_ALIAS.error("testing the logger, error")
    _LOG_ALIAS.critical("testing the logger, info")

    lines = fs.lines_from_file(LOG_FILE)

    assert len(lines) == 3
    assert lines[0].find("WARNING") != -1
    assert lines[1].find("ERROR") != -1
    assert lines[2].find("CRITICAL") != -1

    fs.safe_rm(LOG_FILE)
    assert not os.path.exists(LOG_FILE)


def test_log_data_error_level():
    """Making sure that the logger works at error level."""
    _LOG_ALIAS = lg.init_logger("unit-testing", LOG_FILE, start=True)
    assert not fs.lines_from_file(LOG_FILE)
    _LOG_ALIAS.setLevel(logging.ERROR)
    _LOG_ALIAS.debug("testing the logger, debug")
    _LOG_ALIAS.info("testing the logger, info")
    _LOG_ALIAS.warning("testing the logger, warning")
    _LOG_ALIAS.error("testing the logger, error")
    _LOG_ALIAS.critical("testing the logger, info")

    lines = fs.lines_from_file(LOG_FILE)

    assert len(lines) == 2
    assert lines[0].find("ERROR") != -1
    assert lines[1].find("CRITICAL") != -1

    fs.safe_rm(LOG_FILE)
    assert not os.path.exists(LOG_FILE)


def test_log_data_critical_level():
    """Making sure that the logger works at critical level."""
    _LOG_ALIAS = lg.init_logger("unit-testing", LOG_FILE, start=True)
    assert not fs.lines_from_file(LOG_FILE)
    _LOG_ALIAS.setLevel(logging.CRITICAL)
    _LOG_ALIAS.debug("testing the logger, debug")
    _LOG_ALIAS.info("testing the logger, info")
    _LOG_ALIAS.warning("testing the logger, warning")
    _LOG_ALIAS.error("testing the logger, error")
    _LOG_ALIAS.critical("testing the logger, info")

    lines = fs.lines_from_file(LOG_FILE)

    assert len(lines) == 1
    assert lines[0].find("CRITICAL") != -1

    fs.safe_rm(LOG_FILE)
    assert not os.path.exists(LOG_FILE)
