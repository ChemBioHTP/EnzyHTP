"""Testing enzy_htp._interface.ncaa_library.py.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-10-18
"""
from pathlib import Path

from enzy_htp.structure import Ligand
from enzy_htp._interface.ncaa_library import search_ncaa_parm_file

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/data/"

def test_search_ncaa_parm_file_standard():
    """test search_ncaa_parm_file() using a standard ncaa_lib that all
     file names fit the format."""
    test_ncaa_lib = f"{DATA_DIR}ncaa_lib/"
    test_lig = Ligand(
        residue_idx=1,
        residue_name="H5J",
        atoms=[])

    mol_desc, parm_list = search_ncaa_parm_file(
                            target_res=test_lig,
                            target_method="AM1BCC",
                            ncaa_lib_path=test_ncaa_lib)
    print(mol_desc)
    print(parm_list)

def test_search_ncaa_parm_file_nonstandard():
    """test search_ncaa_parm_file() using a standard ncaa_lib that all
     file names dont fit the format if not necessary."""
    test_ncaa_lib = f"{DATA_DIR}ncaa_lib_non_standard/"
    test_lig = Ligand(
        residue_idx=1,
        residue_name="H5J",
        atoms=[])

    mol_desc, parm_list = search_ncaa_parm_file(
                            target_res=test_lig,
                            target_method="AM1BCC",
                            ncaa_lib_path=test_ncaa_lib)
    print(mol_desc)
    print(parm_list)
