"""Testing enzy_htp._interface.rosetta_interface.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-02-16
"""
from pathlib import Path
import re
import pytest
import os

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp._config.armer_config import ARMerConfig
from enzy_htp import PDBParser
from enzy_htp import interface

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()
ri = interface.rosetta

def test_relax_w_stru():
    """as said in the name. make sure relax works Structure() input."""
    assert False
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    scores = ri.relax(test_stru, 1)
    print(scores)
