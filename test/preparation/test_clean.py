"""Test for enzy_htp.preparation.clean.py for science API related to structure cleaning

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-09-22
"""
import os

import logging
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Structure, PDBParser
from enzy_htp.preparation import clean as cl

# _LOGGER.setLevel(logging.DEBUG)
CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data/" 
WORK_DIR = f"{CURR_DIR}/work_dir/"
sp = PDBParser()

def test_remove_solvent():
    pdb_file_path = f"{DATA_DIR}1Q4T_solvent_cofactor_2chain.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    cl.remove_solvent(stru)
    assert len(stru.residues) == 284

def test_remove_solvent_empty_chain():
    """
    the case there will end up empty chains after remove solvent
    """
    pdb_file_path = f"{DATA_DIR}3NIR.pdb"
    stru: Structure = sp.get_structure(pdb_file_path)
    cl.remove_solvent(stru)
    assert len(stru) == 1


