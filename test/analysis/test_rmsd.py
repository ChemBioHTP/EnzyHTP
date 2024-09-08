#! python3
# -*- encoding: utf-8 -*-
'''
@File    :   test_rmsd.py
@Created :   2024/08/21 20:37
@Author  :   Zhong, Yinjie
@Version :   2.0
@Contact :   yinjie.zhong@vanderbilt.edu
'''

# Here put the import lib.
import glob
import pytest
import os
import numpy as np

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp import PDBParser, interface, _LOGGER

from enzy_htp.analysis import rmsd_of_region, rmsd_with_pattern

from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.core.general import EnablePropagate, get_itself

from enzy_htp._interface.amber_interface import AmberNCParser, AmberRSTParser, AmberMDCRDParser
from enzy_htp.structure.structure_io import PrmtopParser, PDBParser

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()


def test_rmsd_of_region():
    """Test the `rmsd_of_region` function."""
    pdb_file = f"{DATA_DIR}/test_spi.pdb"
    prmtop = f"{DATA_DIR}/test_spi.prmtop"
    data = f"{DATA_DIR}/test_spi.mdcrd"
    parser = PDBParser()

    stru = parser.get_structure(pdb_file)

    mdcrd_parser=AmberMDCRDParser(prmtop).get_coordinates

    structure_ensemble = StructureEnsemble(
        topology=stru,
        top_parser=get_itself,
        coord_parser=mdcrd_parser,
        coordinate_list=data,
    )
    
    rmsd = rmsd_of_region(structure_ensemble=structure_ensemble, ca_only=True)
    _LOGGER.info(f"The RMSD value is {rmsd}.")

    fs.safe_rmdir("scratch")

    assert rmsd < 10

def test_rmsd_with_pattern():
    """Test the `rmsd_with_pattern` function."""
    parser = PDBParser()
    pdb_file = f"{DATA_DIR}/test_spi.pdb"
    prmtop = f"{DATA_DIR}/test_spi.prmtop"
    data = f"{DATA_DIR}/test_spi.mdcrd"

    stru = parser.get_structure(pdb_file)
    ref_structure = parser.get_structure(f"{DATA_DIR}/spi_average_structure.pdb")

    mdcrd_parser=AmberMDCRDParser(prmtop).get_coordinates

    structure_ensemble = StructureEnsemble(
        topology=stru,
        top_parser=get_itself,
        coord_parser=mdcrd_parser,
        coordinate_list=data,
    )
    
    rmsd = rmsd_with_pattern(structure_ensemble=structure_ensemble, reference_structure=ref_structure, mask_pattern="resi 7+11+27+40+41+43+165+168+169+170+199+210+211 not elem H")
    _LOGGER.info(f"The RMSD value is {rmsd}.")

    fs.safe_rmdir("scratch")

    assert rmsd < 10
