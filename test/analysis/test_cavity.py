"""Testing enzy_htp.analysis.cavity

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Created: 2025-06-02
"""

# Here put the import lib.
from os import path
import glob
import pytest
import numpy as np

from enzy_htp.structure import PDBParser
from enzy_htp.analysis.cavity import identify_stru_cavities
import enzy_htp.core.file_system as fs

DATA_DIR = f"{path.dirname(path.abspath(__file__))}/data/"
WORK_DIR = f"{path.dirname(path.abspath(__file__))}/work_dir/"
sp = PDBParser()

def test_identify_stru_cavities():
    """Test `identify_stru_cavities` function."""
    fs.safe_mkdir(WORK_DIR)
    
    pdb_filepath = path.join(DATA_DIR, "cavity_calc", "aclHMT-ETI-SAH_no-ETI.pdb")
    stru = sp.get_structure(pdb_filepath)
    cavities = identify_stru_cavities(stru=stru, work_dir=WORK_DIR)
    fs.safe_rmdir(WORK_DIR)

    assert len(cavities) == 10
    cavity = cavities[0]
    assert len(cavity.boundary_residues) == 4
    assert len(cavity.inner_residues) == 31
    assert abs(cavity.volume - 1147) < 1
    assert abs(cavity.software_report_volume - 947) < 1
