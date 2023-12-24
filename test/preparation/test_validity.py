"""Test for enzy_htp.preparation.validity.py for science API related to structure validity

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-12-19
"""
import os

from enzy_htp.structure import PDBParser
from enzy_htp.preparation import validity as vd

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data/"
WORK_DIR = f"{CURR_DIR}/work_dir/"
sp = PDBParser()

def test_is_structure_valid():
    """test using an example structure"""
    test_pdb = f"{DATA_DIR}1NVG_metalcenter_noligand.pdb"
    test_stru = sp.get_structure(test_pdb)
    vd.is_structure_valid(test_stru, print_report=True)
