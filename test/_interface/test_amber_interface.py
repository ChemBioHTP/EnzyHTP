"""Testing the enzy_htp._interface.AmberInterface class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: QZ Shao <shaoqz@icloud.com>
Date: 2022-06-03
"""
import glob
import os
import shutil
import pytest
from pathlib import Path
from typing import Union

from enzy_htp.core.exception import tLEaPError
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
from enzy_htp.core import file_system as fs
from enzy_htp._interface.amber_interface import (
    AmberParameterizer,
    AmberMDStep)
import enzy_htp.structure as struct
from enzy_htp import interface
from enzy_htp import config as eh_config

MM_BASE_DIR = Path(__file__).absolute().parent
MM_DATA_DIR = f"{MM_BASE_DIR}/data/"
MM_WORK_DIR = f"{MM_BASE_DIR}/work_dir/"
MINIMIZE_INPUT_1 = f"{MM_DATA_DIR}/min_1.inp"
TARGET_MINIMIZE_INPUT_1 = f"{MM_DATA_DIR}/target_min_1.inp"
MINIMIZE_INPUT_2 = f"{MM_DATA_DIR}/min_2.inp"
TARGET_MINIMIZE_INPUT_2 = f"{MM_DATA_DIR}/target_min_2.inp"

ANTECHAMBER_FNAMES = "ANTECHAMBER_AC.AC ANTECHAMBER_AC.AC0 ANTECHAMBER_AM1BCC.AC ANTECHAMBER_AM1BCC_PRE.AC ANTECHAMBER_BOND_TYPE.AC ANTECHAMBER_BOND_TYPE.AC0 ANTECHAMBER_PREP.AC ANTECHAMBER_PREP.AC0 ATOMTYPE.INF NEWPDB.PDB PREP.INF sqm.in sqm.out sqm.pdb".split(
)

# region Tools
def touch(fname: Union[str, Path]):
    if not isinstance(fname, Path):
        fpath = Path(fname)
    else:
        fpath = fname

    if not fpath.exists():
        fh = open(fpath, 'w')
        fh.close()


def exe_exists(exe: str) -> bool:
    """Helper method that checks if an executable exists."""
    return shutil.which(exe) is not None


def files_equivalent(fname1: str, fname2: str) -> bool:
    """Helper method that checks if the contents of two files are equivalent."""
    lines1, lines2 = fs.lines_from_file(fname1), fs.lines_from_file(fname2)

    if len(lines1) != len(lines2):
        return False

    for l1, l2 in zip(lines1, lines2):
        if l1 != l2:
            return False

    return True
# endregion Tools

def test_run_tleap():
    """test function runs as expected"""
    ai = interface.amber
    temp_test_file = f"{MM_BASE_DIR}/work_dir/test_run_tleap.pdb"
    tleap_in_str = f"""source leaprc.protein.ff14SB
a = loadpdb {MM_DATA_DIR}KE_07_R7_2_S.pdb
savepdb a {temp_test_file}
quit
"""

    ai.run_tleap(tleap_in_str)

    assert os.path.exists(temp_test_file)
    with open(temp_test_file) as f:
        assert len(f.readlines()) == 3981
    if _LOGGER.level > 10:
        assert not os.path.exists(f"{eh_config['system.SCRATCH_DIR']}/tleap.out")
        assert not os.path.exists(f"{eh_config['system.SCRATCH_DIR']}/tleap.in")
    fs.safe_rm(temp_test_file)


def test_run_tleap_dev():
    """develop use test function. Dont do any assert"""
    ai = interface.amber
    test_input_pdb = f"{MM_DATA_DIR}KE_07_R7_2_S.pdb"
    temp_test_file = f"{MM_BASE_DIR}/work_dir/test_run_tleap.pdb"
    tleap_in_str = f"""source leaprc.protein.ff14SB
a = loadpdb {test_input_pdb}
savepdb a {temp_test_file}
quit
"""

    ai.run_tleap(tleap_in_str, additional_search_path=["./", "../"])
    fs.safe_rm(temp_test_file)


def test_run_tleap_w_error(caplog):
    """test if the function captures the a errored run of tleap"""
    ai = interface.amber
    temp_test_pdb = f"{MM_BASE_DIR}/work_dir/test_run_tleap.pdb"
    temp_test_prmtop = f"{MM_BASE_DIR}/work_dir/test_run_tleap.prmtop"
    temp_test_inpcrd = f"{MM_BASE_DIR}/work_dir/test_run_tleap.inpcrd"
    # e1
    tleap_in_str = f"""source leaprc.protein.ff24SB
a = loadpdb {MM_DATA_DIR}KE_07_R7_2_S.pdb
savepdb a {temp_test_pdb}
quit
"""
    with pytest.raises(tLEaPError) as e:
        ai.run_tleap(tleap_in_str)
    assert "Could not open file leaprc.protein.ff24SB: not found" in caplog.text
    # e2
    tleap_in_str = f"""source leaprc.protein.ff14SB
a = loadpdb {MM_DATA_DIR}KE_07_R7_2_R.pdb
savepdb a {temp_test_pdb}
quit
"""
    with pytest.raises(tLEaPError) as e:
        ai.run_tleap(tleap_in_str)
    assert "KE_07_R7_2_R.pdb: No such file or directory" in caplog.text
    # e3
    tleap_in_str = f"""source leaprc.protein.ff14SB
a = loadpdb {MM_DATA_DIR}KE_07_R7_2_S.pdb
saveamberparm a {temp_test_prmtop} {temp_test_inpcrd}
quit
"""
    with pytest.raises(tLEaPError) as e:
        ai.run_tleap(tleap_in_str)

    fs.safe_rm(temp_test_pdb)
    fs.safe_rm(temp_test_prmtop)
    fs.safe_rm(temp_test_inpcrd)


def test_find_tleap_error():
    """make sure the function correctly find all the errors in example files"""
    ai = interface.amber
    ERR_EXP_DIR = f"{MM_DATA_DIR}tleap_errors/"
    error_example_and_reason = [
        (f"{ERR_EXP_DIR}e1.out", ["Could not open file leaprc.protein.ff24SB: not found"]),
        (f"{ERR_EXP_DIR}e2.out", ["FATAL:  Atom .R<H5J 254>.A<NAL 2> does not have a type.", "Failed to generate parameters"]),
    ]

    for error_file, error_key_list in error_example_and_reason:
        tleap_error = ai._find_tleap_error(error_file)
        for error_key in error_key_list:
            assert error_key in tleap_error.error_info_str


def test_tleap_clean_up_stru(helpers):
    """make sure the function correctly find all the errors in example files.
    the test file is a truncated KE pdb also deleted the side chain atoms of residue 4
    and changed its name to TRP mimiking a mutation sceniro."""
    ai = interface.amber
    test_input_pdb = f"{MM_DATA_DIR}tleap_clean_up_test_KE.pdb"
    test_out_path = f"{MM_WORK_DIR}tleap_clean_up_out.pdb"
    test_answer_path = f"{MM_DATA_DIR}tleap_clean_up_answer_KE.pdb"

    ai.tleap_clean_up_stru(test_input_pdb, test_out_path, if_align_index=True)

    assert helpers.equiv_files(test_out_path, test_answer_path)
    fs.safe_rm(test_out_path)


def test_build_md_parameterizer_default_value():
    """as said in the name. Assert a charge_method default value
    as a sample."""
    ai = interface.amber
    param_worker: AmberParameterizer = ai.build_md_parameterizer()
    assert param_worker.charge_method == "AM1BCC"


def test_amber_parameterizer_engine():
    """as said in the name."""
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer()
    assert test_param_worker.engine == "Amber"


def test_amber_parameterizer_run_lv_1():
    """level 1 test of the parameterizer.
    Test structure diversity:
    - single polypeptide chain"""
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer(
        ncaa_param_lib_path=f"{MM_DATA_DIR}/ncaa_lib_empty"
    )
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/KE_wo_S.pdb")
    params = test_param_worker.run(test_stru)
    assert params.is_valid()
    for f in params.file_list:
        fs.safe_rm(f)
    fs.safe_rmdir(test_param_worker.parameterizer_temp_dir)
    fs.safe_rmdir(eh_config["system.SCRATCH_DIR"])


def test_amber_parameterizer_run_lv_2():
    """level 2 test of the parameterizer.
    Test structure diversity:
    - single polypeptide chain
    - single substrate (CHON)
    ** use existing parm files for ligand"""
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer(
        ncaa_param_lib_path=f"{MM_DATA_DIR}/ncaa_lib",
        force_fields=[
            "leaprc.protein.ff14SB",
            "leaprc.gaff",
            "leaprc.water.tip3p",
        ]
    )
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb")
    params = test_param_worker.run(test_stru)
    assert params.is_valid()
    for f in params.file_list:
        fs.safe_rm(f)
    fs.safe_rmdir(test_param_worker.parameterizer_temp_dir)
    fs.safe_rmdir(eh_config["system.SCRATCH_DIR"])


def test_amber_parameterizer_run_lv_3():
    """level 3 test of the parameterizer.
    Test structure diversity:
    - single polypeptide chain
    - single substrate (CHON)"""
    temp_mol2_path = f"{MM_DATA_DIR}/ncaa_lib_empty/H5J_AM1BCC-GAFF2.mol2"
    temp_frcmod_path = f"{MM_DATA_DIR}/ncaa_lib_empty/H5J_AM1BCC-GAFF2.frcmod"
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer(
        ncaa_param_lib_path=f"{MM_DATA_DIR}/ncaa_lib_empty"
    )
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    params = test_param_worker.run(test_stru)
    assert params.is_valid()

    for f in params.file_list:
        fs.safe_rm(f)
    fs.safe_rm(temp_mol2_path)
    fs.safe_rm(temp_frcmod_path)
    fs.safe_rmdir(test_param_worker.parameterizer_temp_dir)
    fs.safe_rmdir(eh_config["system.SCRATCH_DIR"])


def test_amber_parameterizer_run_lv_4():
    """level 4 test of the parameterizer.
    Test structure diversity:
    - 2 polypeptide chain
    - 2 substrate (CHN, CHONP), special net charge
    * use lib to make test faster"""
    temp_mol2_path = f"{MM_DATA_DIR}/ncaa_lib_empty/PUT_AM1BCC-GAFF.mol2"
    temp_frcmod_path = f"{MM_DATA_DIR}/ncaa_lib_empty/PUT_AM1BCC-GAFF.frcmod"
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer(
        ncaa_param_lib_path=f"{MM_DATA_DIR}/ncaa_lib_empty",
        force_fields=[
            "leaprc.protein.ff14SB",
            "leaprc.gaff",
            "leaprc.water.tip3p",
        ]
    )
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/puo_put.pdb")
    test_stru.assign_ncaa_chargespin({"PUT" : (1,1),
                                      "FAD" : (0,1)})
    params = test_param_worker.run(test_stru)
    assert params.is_valid()

    for f in params.file_list:
        fs.safe_rm(f)
    fs.safe_rm(temp_mol2_path)
    fs.safe_rm(temp_frcmod_path)
    fs.safe_rmdir(test_param_worker.parameterizer_temp_dir)
    fs.safe_rmdir(eh_config["system.SCRATCH_DIR"])


def test_amber_parameterizer_run_lv_5(): #TODO
    """level 5 test of the parameterizer.
    Test structure diversity:
    - 2 polypeptide chain
    - 1 substrate (CHONP)
    - 1 modified amino acid (CHONP)"""
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer(
        ncaa_param_lib_path=f"{MM_DATA_DIR}/ncaa_lib_empty"
    )
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/3cfr-slp-pea_ah.pdb")
    test_param_worker.run(test_stru)


def test_amber_parameterizer_run_lv_6(): #TODO
    """level 6 test of the parameterizer.
    Test structure diversity:
    - 2 polypeptide chain
    - 1 substrate (CHONP)
    - 1 modified amino acid (CHONP)"""
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer(
        ncaa_param_lib_path=f"{MM_DATA_DIR}/ncaa_lib_empty"
    )
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/tyna_clean.pdb")
    test_param_worker.run(test_stru)


def test_check_gaff_type():
    """test _check_gaff_type in AmberParameterizer"""
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer(
        ncaa_param_lib_path=f"{MM_DATA_DIR}/ncaa_lib_empty",
        force_fields=[
            "leaprc.protein.ff14SB",
            "leaprc.gaff2",
            "leaprc.water.tip3p",
        ]
    )
    assert test_param_worker._check_gaff_type() == "GAFF2"

    test_param_worker.force_fields = [
            "leaprc.protein.ff14SB",
            "leaprc.gaff",
            "leaprc.water.tip3p",
        ]
    assert test_param_worker._check_gaff_type() == "GAFF"

    test_param_worker.force_fields = [
            "leaprc.protein.ff14SB",
            "leaprc.water.tip3p",
        ]
    assert test_param_worker._check_gaff_type() is None


def test_run_parmchk2():
    """test the function works well"""
    ai = interface.amber
    test_in_file = f"{MM_DATA_DIR}/ncaa_lib/H5J_AM1BCC-GAFF.prepin"
    temp_frcmod_file = f"{MM_WORK_DIR}/frcmod"

    ai.run_parmchk2(test_in_file, temp_frcmod_file, gaff_type="GAFF")

    assert os.path.exists(temp_frcmod_file)
    assert len(fs.lines_from_file(temp_frcmod_file)) == 23
    fs.safe_rm(temp_frcmod_file)


def test_run_antechamber():
    """test the function works well"""
    ai = interface.amber
    test_in_file = f"{MM_DATA_DIR}/H5J.pdb"
    temp_mol2_file = f"{MM_WORK_DIR}/H5J.mol2"

    ai.run_antechamber(test_in_file, temp_mol2_file,
                        net_charge=0, spin=1,
                        charge_method="bcc")

    assert os.path.exists(temp_mol2_file)
    assert len(fs.lines_from_file(temp_mol2_file)) == 43
    for temp_files in ["ATOMTYPE.INF","NEWPDB.PDB","PREP.INF","sqm.pdb","sqm.in","sqm.out"]:
        assert not os.path.exists(temp_files)
    assert not glob.glob("ANTECHAMBER*")
    fs.safe_rm(temp_mol2_file)


def test_run_antechamber_wrong(caplog):
    """test the function works well reporting the error"""
    ai = interface.amber
    test_in_file = f"{MM_DATA_DIR}/H5J.pdb"
    temp_mol2_file = f"{MM_WORK_DIR}/H5J.mol2"

    with EnablePropagate(_LOGGER):
        with pytest.raises(ValueError) as exe:
            ai.run_antechamber(test_in_file, temp_mol2_file,
                                net_charge=0, spin=1,
                                charge_method="cm1")
        assert "found unsupported charge method" in caplog.text

        with pytest.raises(ValueError) as exe:
            ai.run_antechamber(test_in_file, temp_mol2_file,
                                net_charge=0, spin=1,
                                charge_method="rc")
        assert "no charge_file is provided" in caplog.text

        with pytest.raises(ValueError) as exe:
            ai.run_antechamber(test_in_file, temp_mol2_file,
                                net_charge=0, spin=1,
                                charge_method="resp")
        assert "not gesp or gout" in caplog.text


def test_antechamber_ncaa_to_moldesc():
    """test the function works well"""
    temp_mol2_path = f"{MM_WORK_DIR}/H5J_AM1BCC.mol2"
    ai = interface.amber
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb").ligands[0]
    test_stru.net_charge = 0
    test_stru.spin = 1

    ai.antechamber_ncaa_to_moldesc(ncaa=test_stru,
                                   out_path=temp_mol2_path,
                                   gaff_type="GAFF",
                                   charge_method="AM1BCC",)

    assert os.path.exists(temp_mol2_path)
    assert len(fs.lines_from_file(temp_mol2_path)) == 43
    assert not os.path.exists("scratch/H5J.pdb")
    fs.safe_rm(temp_mol2_path)


def test_antechamber_ncaa_to_moldesc_wrong(caplog):
    """test the function works well"""
    temp_mol2_path = f"{MM_WORK_DIR}/H5J_AM1BCC.mol2"
    ai = interface.amber
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb").ligands[0]

    with EnablePropagate(_LOGGER):
        with pytest.raises(ValueError) as exe:
            ai.antechamber_ncaa_to_moldesc(ncaa=test_stru,
                                        out_path=temp_mol2_path,
                                        gaff_type="GAFF",
                                        charge_method="AM1BCC",)
            assert "does not have charge and spin." in caplog.text


def test_antechamber_ncaa_to_moldesc_resp(): #TODO
    """test the function works well"""
    temp_mol2_path = f"{MM_WORK_DIR}/H5J_RESP.mol2"
    ai = interface.amber
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb").ligands[0]
    test_stru.net_charge = 0
    test_stru.spin = 1

    ai.antechamber_ncaa_to_moldesc(ncaa=test_stru,
                                   out_path=temp_mol2_path,
                                   gaff_type="GAFF",
                                   charge_method="RESP",)

    assert os.path.exists(temp_mol2_path)
    assert len(fs.lines_from_file(temp_mol2_path)) == 43
    assert not os.path.exists("scratch/H5J.pdb")
    fs.safe_rm(temp_mol2_path)


def test_build_md_step_default():
    """as said in the name. Assert a several default value
    as a sample."""
    ai = interface.amber
    md_step: AmberMDStep = ai.build_md_step(length=0.1, core_type="cpu")
    assert md_step.temperature == 300.0
    assert md_step.core_type == "cpu"
    assert md_step.cluster_job_config["res_setting"]["core_type"] == "cpu"
    assert md_step.length == 0.1
    assert md_step.record_period == 0.0001

# region TODO

def test_parse_fmt_uppercase():
    """Testing that the AmberInterface.parse_fmt() method works correctly for uppercase inputs."""
    ai = interface.amber
    assert ai.parse_fmt('%FORMAT(20A4)') == (str, 4)

    assert ai.parse_fmt('%FORMAT(5E16.8)') == (float, 16)

    assert ai.parse_fmt('%FORMAT(10I8)') == (int, 8)


def test_parse_fmt_lowercase():
    """Testing that the AmberInterface.parse_fmt() method works correctly for lowercase inputs."""
    ai = interface.amber
    assert ai.parse_fmt('%FORMAT(20a4)') == (str, 4)

    assert ai.parse_fmt('%FORMAT(5e16.8)') == (float, 16)

    assert ai.parse_fmt('%FORMAT(10i8)') == (int, 8)


def test_parse_fmt_bad_input():
    """Testing that the AmberInterface.parse_fmt() method returns the correct (None,-1) for bad inputs."""

    ai = interface.amber

    assert ai.parse_fmt('') == (None, -1)

    assert ai.parse_fmt('20a4') == (None, -1)

    assert ai.parse_fmt('FORMAT(20a4)') == (None, -1)

    assert ai.parse_fmt('%FORMAT(20a4') == (None, -1)



def test_parse_prmtop():
    """Testing that AmberInterface.parse_prmtop() method works. It essentially loads the supplied
    prmtop file and ensures that the correct number of items are extracted for the majority of keys.
    """
    ai = interface.amber
    data = ai.parse_prmtop(f"{MM_DATA_DIR}/prmtop_1")

    assert data

    n_atoms: int = 23231
    num_bond: int = 85
    num_ang: int = 186
    n_ptra: int = 190
    n_types: int = 19
    n_bond_h: int = 21805
    n_bond_a: int = 1450
    n_thet: int = 3319
    n_phi: int = 6138
    n_phi_h: int = 6551
    nnb: int = 43036
    n_phb: int = 1
    n_res: int = 6967

    assert len(data['ATOM_NAME']) == n_atoms
    assert len(data['CHARGE']) == n_atoms
    assert len(data['ATOMIC_NUMBER']) == n_atoms
    assert len(data['MASS']) == n_atoms
    assert len(data['ATOM_TYPE_INDEX']) == n_atoms
    assert len(data['NUMBER_EXCLUDED_ATOMS']) == n_atoms
    assert len(data['NONBONDED_PARM_INDEX']) == n_types**2
    assert len(data['RESIDUE_LABEL']) == n_res
    assert len(data['RESIDUE_POINTER']) == n_res
    assert len(data['BOND_FORCE_CONSTANT']) == num_bond
    assert len(data['BOND_EQUIL_VALUE']) == num_bond
    assert len(data['ANGLE_FORCE_CONSTANT']) == num_ang
    assert len(data['ANGLE_EQUIL_VALUE']) == num_ang
    assert len(data['DIHEDRAL_FORCE_CONSTANT']) == n_ptra
    assert len(data['DIHEDRAL_PERIODICITY']) == n_ptra
    assert len(data['DIHEDRAL_PHASE']) == n_ptra
    assert len(data['SCEE_SCALE_FACTOR']) == n_ptra
    assert len(data['SOLTY']) == 52
    assert len(data['LENNARD_JONES_ACOEF']) == (n_types * (n_types + 1)) / 2
    assert len(data['LENNARD_JONES_BCOEF']) == (n_types * (n_types + 1)) / 2
    assert len(data['BONDS_INC_HYDROGEN']) == 3 * n_bond_h
    assert len(data['BONDS_WITHOUT_HYDROGEN']) == 3 * n_bond_a
    assert len(data['ANGLES_INC_HYDROGEN']) == 4 * n_thet
    assert len(data['DIHEDRALS_INC_HYDROGEN']) == 5 * n_phi_h
    assert len(data['DIHEDRALS_WITHOUT_HYDROGEN']) == 5 * n_phi
    assert len(data['EXCLUDED_ATOMS_LIST']) == nnb
    assert len(data['HBOND_ACOEF']) == n_phb
    assert len(data['HBOND_BCOEF']) == n_phb
    assert len(data['HBCUT']) == n_phb
    assert len(data['AMBER_ATOM_TYPE']) == n_atoms
    assert len(data['TREE_CHAIN_CLASSIFICATION']) == n_atoms
    assert len(data['JOIN_ARRAY']) == n_atoms
    assert len(data['IROTAT']) == n_atoms
    assert len(data['RADII']) == n_atoms


def test_parse_prmtop_no_file():
    """Checking that AmberInteface.parse_prmtop() throws an error if the listed file does not exist."""

    dne = Path('./dne')

    assert not dne.exists()

    ai = interface.amber
    with pytest.raises(SystemExit) as exe:
        ai.parse_prmtop(dne)

    assert exe


def test_add_charges():
    """Checking that the AmberInterface.add_charges() method works correctly."""

    ai = interface.amber
    ss = struct.PDBParser.get_structure(f"{MM_BASE_DIR}/data/1h1d.pdb")
    assert not ss.has_charges()
    ai.add_charges(ss, f"{MM_BASE_DIR}/data/1h1d_prmtop")
    assert ss.has_charges()


def test_add_charges_bad_file():
    """Checking that the AmberInterface.add_charge() method throws an error for bad inputs."""

    ai = interface.amber
    ss = struct.PDBParser.get_structure(f"{MM_BASE_DIR}/data/1h1d.pdb")

    with pytest.raises(SystemExit) as exe:
        ai.add_charges(ss, f"{MM_BASE_DIR}/data/bad_label_prmtop")

    assert exe

# endregion TODO
