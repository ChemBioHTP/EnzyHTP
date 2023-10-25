"""Testing enzy_htp._interface.ncaa_library.py.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-10-18
"""
import logging
from pathlib import Path

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Ligand
from enzy_htp._interface.ncaa_library import search_ncaa_parm_file

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/data/"
_LOGGER.propagate = True

def test_search_ncaa_parm_file_standard():
    """test search_ncaa_parm_file() using a standard ncaa_lib that all
     file names fit the format. use cache is turned off for the test."""
    test_ncaa_lib = f"{DATA_DIR}ncaa_lib/"
    test_lig = Ligand(
        residue_idx=1,
        residue_name="H5J",
        atoms=[])

    mol_desc, parm_list = search_ncaa_parm_file(
                            target_res=test_lig,
                            target_method="AM1BCC",
                            ncaa_lib_path=test_ncaa_lib,
                            use_cache=False)
    assert fs.base_file_name(mol_desc) == "H5J_AM1BCC"
    assert fs.get_file_ext(mol_desc) == ".prepin"
    assert len(parm_list) == 2

def test_search_ncaa_parm_file_nonstandard():
    """test search_ncaa_parm_file() using a nonstandard ncaa_lib that all
     file names dont fit the format if not necessary."""
    test_ncaa_lib = f"{DATA_DIR}ncaa_lib_non_standard/"
    test_lig = Ligand(
        residue_idx=1,
        residue_name="H5J",
        atoms=[])

    mol_desc, parm_list = search_ncaa_parm_file(
                            target_res=test_lig,
                            target_method="AM1BCC",
                            ncaa_lib_path=test_ncaa_lib,
                            use_cache=False)
    assert fs.base_file_name(mol_desc) == "RES_AM1BCC"
    assert fs.get_file_ext(mol_desc) == ".prepin"
    assert not parm_list

def test_search_ncaa_parm_file_wrongname(caplog):
    """test search_ncaa_parm_file() using a ncaa_lib with wrong filename"""
    test_ncaa_lib = f"{DATA_DIR}ncaa_lib_wrong/"
    test_lig = Ligand(
        residue_idx=1,
        residue_name="H5J",
        atoms=[])

    mol_desc, parm_list = search_ncaa_parm_file(
                            target_res=test_lig,
                            target_method="AM1BCC",
                            ncaa_lib_path=test_ncaa_lib,
                            use_cache=False)
    assert not mol_desc
    assert not parm_list
    assert "NCAA lib file 'method' (H5J) does not match any supported ones." in caplog.text

def test_search_ncaa_parm_file_any():
    """test search_ncaa_parm_file() using a ncaa_lib with wrong filename
    but matched using any and resname from the file"""
    test_ncaa_lib = f"{DATA_DIR}ncaa_lib_wrong/"
    test_lig = Ligand(
        residue_idx=1,
        residue_name="H5J",
        atoms=[])

    mol_desc, parm_list = search_ncaa_parm_file(
                            target_res=test_lig,
                            target_method="any",
                            ncaa_lib_path=test_ncaa_lib,
                            use_cache=False)
    assert mol_desc
    assert not parm_list

def test_search_ncaa_parm_file_cache(caplog):
    """test search_ncaa_parm_file() using a standard ncaa_lib with cache loaded"""
    test_ncaa_lib = f"{DATA_DIR}ncaa_lib/"
    test_lig = Ligand(
        residue_idx=1,
        residue_name="H5J",
        atoms=[])

    search_ncaa_parm_file(target_res=test_lig,
                          target_method="any",
                          ncaa_lib_path=test_ncaa_lib,
                          use_cache=True)
    _LOGGER.setLevel(logging.DEBUG)
    mol_desc, parm_list = search_ncaa_parm_file(
                            target_res=test_lig,
                            target_method="AM1BCC",
                            ncaa_lib_path=test_ncaa_lib,
                            use_cache=True)
    assert "using cache" in caplog.text

    assert fs.base_file_name(mol_desc) == "H5J_AM1BCC"
    assert fs.get_file_ext(mol_desc) == ".prepin"
    assert len(parm_list) == 2
