"""Testing the enzy_htp.structure.structure_operation.general.py

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-22
"""

import os
from enzy_htp.structure import PDBParser, Structure
import enzy_htp.structure.structure_operation as stru_oper

#_LOGGER.setLevel(logging.DEBUG)
CURRDIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f'{CURRDIR}/../data/'
sp = PDBParser()

def test_remove_solvent():
    pdb_file_path = f'{DATA_DIR}1Q4T_ligand_test.pdb'
    stru: Structure = sp.get_structure(pdb_file_path)

    stru_oper.remove_solvent(stru)
    assert len(stru.solvents) == 0

def test_remove_empty_chain():
    pdb_file_path = f'{DATA_DIR}1Q4T_ligand_test.pdb'
    stru: Structure = sp.get_structure(pdb_file_path)

    stru[1].children = []
    stru_oper.remove_empty_chain(stru)
    assert len(stru) == 3

