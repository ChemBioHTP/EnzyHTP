import pytest
import os
from enzy_htp import interface
from enzy_htp.structure import Structure, StructureEnsemble
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp.analysis import residue_pka

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
INT_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../_interface/data/"
STRU_DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../test_data/diversed_stru/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"

def test_residue_pka_structure(pdb_path):
    pdb_file = f"{DATA_DIR}/test_spi.pdb"
    stru = PDBParser().get_structure(pdb_file)


def test_residue_pka_structure_ensemble(pdb_path):
    stru_esm = interface.amber.load_traj(
        prmtop_path=f"{INT_DATA_DIR}/mmpbsa_test_sol.prmtop",
        traj_path=f"{INT_DATA_DIR}/mmpbsa_test_sol_10f.nc",
        ref_pdb=f"{INT_DATA_DIR}/mmpbsa_test_sol.pdb"
    )
