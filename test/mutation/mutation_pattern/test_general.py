"""Testing the enzy_htp.mutation.mutation_pattern.general.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-01-26
"""

import pytest
import os

import enzy_htp.mutation.mutation_pattern.general as m_p
from enzy_htp import PDBParser

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()

@pytest.mark.interface
def test_decode_mutation_pattern():
    """dev run of the function"""
    test_mutation_pattern = "RA154W, r:100[resi 289 around 4 and not resi 36:larger, proj(id 1000, id 2023, positive, 10):more_negative_charge]"
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)

    m_p.decode_mutation_pattern(test_stru, test_mutation_pattern)

def test_seperate_sections():
    """test the function use a made up mutation_pattern for KE
    use the len of the seperation result as the fingerprint to assert"""
    test_mutation_pattern = "RA154W, r:100[resi 289 around 4 and not resi 36:larger, proj(id 1000, id 2023, positive, 10):more_negative_charge], a:[resi 1:all]"
    assert len(m_p.seperate_sections(test_mutation_pattern)) == 3

