"""Testing the enzy_htp._interface.AmberInterface class.
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Author: QZ Shao <shaoqz@icloud.com>
Date: 2022-06-03
"""
import glob
import os
import re
import shutil
from subprocess import CompletedProcess
import pytest
import numpy as np
from pathlib import Path
from typing import Union

from enzy_htp.core.clusters.accre import Accre
from enzy_htp.core.exception import tLEaPError, AmberMDError
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.general import EnablePropagate
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig
from enzy_htp.core import file_system as fs
from enzy_htp._interface.amber_interface import (
    AmberParameterizer,
    AmberParameter,
    AmberMDStep,
    AmberMDResultEgg,)
import enzy_htp.structure as struct
from enzy_htp.structure.structure_constraint import (
    StructureConstraint,
    create_cartesian_freeze,
    create_backbone_freeze,
    create_distance_constraint,
    create_angle_constraint,)
from enzy_htp import interface
from enzy_htp import config as eh_config
from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.structure.structure_selection.general import select_stru

MM_BASE_DIR = Path(__file__).absolute().parent
MM_DATA_DIR = f"{MM_BASE_DIR}/data/"
STRU_DATA_DIR = f"{MM_BASE_DIR}/../structure/data/"
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
    with EnablePropagate(_LOGGER):
        with pytest.raises(tLEaPError) as e:
            ai.run_tleap(tleap_in_str)
        assert "Could not open file leaprc.protein.ff24SB: not found" in caplog.text
    # e2
    tleap_in_str = f"""source leaprc.protein.ff14SB
a = loadpdb {MM_DATA_DIR}KE_07_R7_2_R.pdb
savepdb a {temp_test_pdb}
quit
"""
    with EnablePropagate(_LOGGER):
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


def test_amber_parameterizer_run_lv_1_add_lines():
    """level 1 test of the parameterizer with additional custom lines to tleap.in.
    Test structure diversity:
    - single polypeptide chain"""
    ai = interface.amber
    test_param_worker: AmberParameterizer = ai.build_md_parameterizer(
        ncaa_param_lib_path=f"{MM_DATA_DIR}/ncaa_lib_empty",
        additional_tleap_lines=["#this is a test line"],
        keep_tleap_in=True,
    )
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/KE_wo_S.pdb")
    params = test_param_worker.run(test_stru)

    # NOTE potential failed test due to a different temp path is generated
    with open(f"{eh_config['system.SCRATCH_DIR']}/tleap.in") as f:
        assert "#this is a test line" in f.read()

    for f in params.file_list:
        fs.safe_rm(f)
    fs.safe_rmdir(test_param_worker.parameterizer_temp_dir)
    fs.safe_rmdir(eh_config["system.SCRATCH_DIR"])


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
    test_stru.multiplicity = 1

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
    assert Exception("TODO")
    temp_mol2_path = f"{MM_WORK_DIR}/H5J_RESP.mol2"
    ai = interface.amber
    test_stru = struct.PDBParser().get_structure(
        f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb").ligands[0]
    test_stru.net_charge = 0
    test_stru.multiplicity = 1

    ai.antechamber_ncaa_to_moldesc(ncaa=test_stru,
                                   out_path=temp_mol2_path,
                                   gaff_type="GAFF",
                                   charge_method="RESP",)

    assert os.path.exists(temp_mol2_path)
    assert len(fs.lines_from_file(temp_mol2_path)) == 43
    assert not os.path.exists("scratch/H5J.pdb")
    fs.safe_rm(temp_mol2_path)


def test_build_md_step_default():
    """as said in the name. Assert several default values
    as samples."""
    ai = interface.amber
    md_step: AmberMDStep = ai.build_md_step(length=0.1, core_type="cpu")
    assert md_step.temperature == 300.0
    assert md_step.core_type == "cpu"
    assert md_step.cluster_job_config["res_keywords"]["core_type"] == "cpu"
    assert md_step.length == 0.1
    assert md_step.record_period == 0.0001


def test_build_md_step_res_keywords():
    """as said in the name. Assert several default values
    as samples."""
    ai = interface.amber
    md_step: AmberMDStep = ai.build_md_step(length=0.1, core_type="cpu",
                                            cluster_job_config={
                                                "cluster" : None,
                                                "res_keywords" : {"partition" : "production",
                                                                  "account" : "yang_lab",}
                                            })
    assert md_step.cluster_job_config["cluster"] is None
    assert md_step.cluster_job_config["res_keywords"]["core_type"] == "cpu"
    assert md_step.cluster_job_config["res_keywords"]["partition"] == "production"
    assert md_step.cluster_job_config["res_keywords"]["nodes"] == "1"
    assert md_step.cluster_job_config["res_keywords"]["node_cores"] ==  "16"
    assert md_step.cluster_job_config["res_keywords"]["job_name"] ==  "MD_EnzyHTP"
    assert md_step.cluster_job_config["res_keywords"]["mem_per_core"] ==  "3G"
    assert md_step.cluster_job_config["res_keywords"]["walltime"] ==  "1-00:00:00"
    assert md_step.cluster_job_config["res_keywords"]["account"] ==  "yang_lab"


def test_write_to_mdin_from_raw_dict():
    """test to make sure _write_to_mdin_from_raw_dict() works as expected.
    using a dict from old EnzyHTP Class_Conf.Amber.conf_heat as an example"""
    test_raw_dict = {
        'title': 'Heat',
        'namelists': [
           {'type': 'cntrl',
            'config': {
                'imin': 0, 'ntx': 1, 'irest': 0,
                'ntc': 2, 'ntf': 2,
                'cut': 10.0,
                'nstlim': 20000, 'dt': 0.002,
                'tempi': 0.0, 'temp0': 300.0,
                'ntpr': 200, 'ntwx': 20000,
                'ntt': 3, 'gamma_ln': 5.0,
                'ntb': 1, 'ntp':0,
                'iwrap': 1,
                'nmropt': 1,
                'ig': -1,
                'ntr': 1, 'restraint_wt': 2.0, 'restraintmask': "'@C,CA,N'",
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'TEMP0'",
                'istep1': 0, 'istep2': 18000,
                'value1': 0.0, 'value2': 300.0,
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'TEMP0'",
                'istep1': 18001, 'istep2': 20000,
                'value1': 300.0, 'value2': 300.0,
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'END'",
                }
            },
        ],
        'file_redirection': {
            'DISANG': {
                'path': './MD//0.rs',
                'content': [
                    {'ialtd': 0,
                    'iat': [3976, 1579],
                    'r1': 2.150,
                    'r2': 2.350,
                    'r3': 2.450,
                    'r4': 2.650,
                    'rk2': 200.0,
                    'rk3': 200.0},
                    {'ialtd': 0,
                    'iat': [3975,3976,1579],
                    'r1': 150.0,
                    'r2': 170.0,
                    'r3': 190.0,
                    'r4': 210.0,
                    'rk2': 200.0,
                    'rk3': 200.0}
                ]}
        },
        'group_info': [],
    }
    test_temp_mdin = f"{MM_WORK_DIR}/test_mdin_from_raw_dict.in"
    answer_temp_mdin = f"{MM_DATA_DIR}/answer_mdin_from_raw_dict.in"
    fs.safe_mkdir("MD/")
    ai = interface.amber
    ai._write_to_mdin_from_raw_dict(test_raw_dict, test_temp_mdin)
    assert files_equivalent(test_temp_mdin, answer_temp_mdin)
    fs.safe_rm(test_temp_mdin)
    fs.safe_rm("MD/0.rs")
    fs.safe_rmdir("MD/")


def test_write_disang_file():
    """test using an example raw dict list"""
    test_dict_list = [
        {'ialtd': 0,
        'iat': [3976, 1579],
        'r1': 2.150,
        'r2': 2.350,
        'r3': 2.450,
        'r4': 2.650,
        'rk2': 200.0,
        'rk3': 200.0},
        {'ialtd': 0,
        'iat': [3975,3976,1579],
        'r1': 150.0,
        'r2': 170.0,
        'r3': 190.0,
        'r4': 210.0,
        'rk2': 200.0,
        'rk3': 200.0}]
    test_disang = f"{MM_WORK_DIR}/test_disang_from_raw_dict.rs"
    answer_disang = f"{MM_DATA_DIR}/answer_disang_from_raw_dict.rs"
    ai = interface.amber
    ai.write_disang_file(test_dict_list, test_disang)
    assert files_equivalent(test_disang, answer_disang)

    fs.safe_rm(test_disang)


def test_parse_md_config_dict_to_raw_wo_cons():
    """test to make sure _parse_md_config_dict_to_raw() works as expected.
    using a dict from old EnzyHTP Class_Conf.Amber.conf_heat as an example"""
    answer_raw_dict = {
        'title': 'Heat',
        'namelists': [
           {'type': 'cntrl',
            'config': {
                'imin': 0, 'ntx': 1, 'irest': 0,
                'ntc': 2, 'ntf': 2,
                'cut': 10.0,
                'nstlim': 20000, 'dt': 0.002,
                'tempi': 0.0, 'temp0': 300.0,
                'ntpr': 200, 'ntwx': 200,
                'ntt': 3, 'gamma_ln': 5.0,
                'ntb': 1, 'ntp':0,
                'iwrap': 1,
                'ig': -1,
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'TEMP0'",
                'istep1': 0, 'istep2': 18000,
                'value1': 0.0, 'value2': 300.0,
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'TEMP0'",
                'istep1': 18001, 'istep2': 20000,
                'value1': 300.0, 'value2': 300.0,
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'END'",
                }
            },
        ],
        'file_redirection': {},
        'group_info': [],
    }
    test_md_config_dict = {
            "name" : "Heat",
            "length" : 0.04, # ns
            "timestep" : 0.000002, # ns
            "minimize" : False,
            "temperature" : [(0.0, 0.0), (0.036, 300.0), (0.04, 300.0)],
            "thermostat" : "langevin",
            "pressure_scaling" : "none",
            "constrain" : None,
            "restart" : False,
            "if_report" : True,
            "record_period" : 0.0004, # ns
            "mdstep_dir" : "./MD",
    }

    ai = interface.amber
    test_raw_dict = ai._parse_md_config_dict_to_raw(test_md_config_dict)
    assert test_raw_dict == answer_raw_dict


def test_parse_md_config_dict_to_raw_minimize():
    """test to make sure _parse_md_config_dict_to_raw() works as expected.
    using a dict from old EnzyHTP Class_Conf.Amber.conf_min as an example"""
    answer_raw_dict = {
        'title': 'Min',
        'namelists': [
           {'type': 'cntrl',
            'config': {
                'imin': 1, 'ntx': 1, 'irest': 0,
                'ntc': 2, 'ntf': 2,
                'cut': 10.0,
                'maxcyc': 20000, 'ncyc': 10000,
                'ntpr': 200, 'ntwx': 0,
                }
            },
        ],
        'file_redirection': {},
        'group_info': [],
    }
    test_md_config_dict = {
            "name" : "Min",
            "length" : 20000, # cycle
            "timestep" : 0.000002, # ns
            "minimize" : True,
            "temperature" : 300.0,
            "thermostat" : "langevin",
            "pressure_scaling" : "none",
            "constrain" : None,
            "restart" : False,
            "if_report" : True,
            "record_period" : 0.0004, # ns
            "mdstep_dir" : "./MD",
    }

    ai = interface.amber
    test_raw_dict = ai._parse_md_config_dict_to_raw(test_md_config_dict)
    assert test_raw_dict == answer_raw_dict


def test_parse_md_config_dict_to_raw_w_cons():
    """test to make sure _parse_md_config_dict_to_raw() works as expected.
    using a dict from old EnzyHTP Class_Conf.Amber.conf_heat as an example"""
    test_pdb = f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb"
    test_stru = struct.PDBParser().get_structure(test_pdb)
    answer_raw_dict = {
        'title': 'Heat',
        'namelists': [
           {'type': 'cntrl',
            'config': {
                'imin': 0, 'ntx': 1, 'irest': 0,
                'ntc': 2, 'ntf': 2,
                'cut': 10.0,
                'nstlim': 20000, 'dt': 0.002,
                'tempi': 0.0, 'temp0': 300.0,
                'ntpr': 200, 'ntwx': 200,
                'ntt': 3, 'gamma_ln': 5.0,
                'ntb': 1, 'ntp':0,
                'iwrap': 1,
                'nmropt': 1,
                'ig': -1,
                'ntr': 1, 'restraint_wt': 2.0, 'restraintmask': "'@C,CA,N'",
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'TEMP0'",
                'istep1': 0, 'istep2': 18000,
                'value1': 0.0, 'value2': 300.0,
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'TEMP0'",
                'istep1': 18001, 'istep2': 20000,
                'value1': 300.0, 'value2': 300.0,
                }
            },
           {'type': 'wt',
            'config': {
                'type': "'END'",
                }
            },
        ],
        'file_redirection': {
            'DISANG': {
                'path': './MD//0.rs',
                'content': [
                    {'ialtd': 0,
                    'iat': [3976, 1579],
                    'r1': 2.150,
                    'r2': 2.350,
                    'r3': 2.450,
                    'r4': 2.650,
                    'rk2': 200.0,
                    'rk3': 200.0},
                    {'ialtd': 0,
                    'iat': [3975,3976,1579],
                    'r1': 150.0,
                    'r2': 170.0,
                    'r3': 190.0,
                    'r4': 210.0,
                    'rk2': 200.0,
                    'rk3': 200.0}
                ]} 
        },
        'group_info': [],
    }
    test_md_config_dict = {
        "name" : "Heat",
        "minimize" : False,
        "pressure_scaling" : "none",
        "restart" : False,
        "length" : 0.04, # ns
        "timestep" : 0.000002, # ns
        "temperature" : [(0.0, 0.0), (0.036, 300.0), (0.04, 300.0)],
        "thermostat" : "langevin",
        "constrain" : [
            create_backbone_freeze(test_stru),
            create_distance_constraint(
                "B.254.H2", "A.101.OE2", 2.4, test_stru),
            create_angle_constraint(
                "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru)],
        "if_report" : True,
        "record_period" : 0.0004,
        "mdstep_dir" : "./MD",
    }

    ai = interface.amber
    test_raw_dict = ai._parse_md_config_dict_to_raw(test_md_config_dict)
    assert test_raw_dict == answer_raw_dict


def test_amber_md_step_make_job():
    """test to make sure AmberMDStep.make_job() works as expected.
    w/o constraint."""
    ai = interface.amber
    md_step = ai.build_md_step(length=0.1) # 300K, NPT by default
    test_inpcrd = f"{MM_DATA_DIR}/KE_07_R7_S.inpcrd"
    test_prmtop = f"{MM_DATA_DIR}/KE_07_R7_S.prmtop"
    test_params = AmberParameter(test_inpcrd, test_prmtop)
    test_job, test_md_egg = md_step.make_job(test_params)
    
    answer_pattern = r"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=MD_EnzyHTP
#SBATCH --partition=<fillthis>
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --account=<fillthis>
#SBATCH --export=NONE
(?:#SBATCH --exclude=.+)?

# Script generated by EnzyHTP [0-9]\.[0-9]\.[0-9] in [0-9]+-[0-9]+-[0-9]+ [0-9]+:[0-9]+:[0-9]+

source /home/shaoq1/bin/amber_env/amber22\.sh

pmemd\.cuda -O -i \./MD/amber_md_step_?[0-9]*\.in -o \./MD/amber_md_step\.out -p .*test/_interface/data//KE_07_R7_S\.prmtop -c .*test/_interface/data//KE_07_R7_S\.inpcrd -r \./MD/amber_md_step\.rst -ref .*test/_interface/data//KE_07_R7_S\.inpcrd -x \./MD/amber_md_step\.nc 
"""
    assert re.match(answer_pattern, test_job.sub_script_str)
    assert test_md_egg.traj_path == './MD/amber_md_step.nc'
    assert test_md_egg.traj_log_path == './MD/amber_md_step.out'
    assert test_md_egg.rst_path == './MD/amber_md_step.rst'
    assert Path(test_md_egg.prmtop_path) == Path('test/_interface/data//KE_07_R7_S.prmtop').absolute()
    fs.safe_rmdir(md_step.work_dir)


def test_amber_md_step_make_job_w_cons():
    """test to make sure AmberMDStep.make_job() works as expected.
    w/ constraint. Just make sure no exceptions are raised"""
    test_inpcrd = f"{MM_DATA_DIR}/KE_07_R7_S.inpcrd"
    test_prmtop = f"{MM_DATA_DIR}/KE_07_R7_S.prmtop"
    test_pdb = f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb"
    test_stru = struct.PDBParser().get_structure(test_pdb)
    test_constrains = [
        create_backbone_freeze(test_stru),
        create_distance_constraint(
            "B.254.H2", "A.101.OE2", 2.4, test_stru),
        create_angle_constraint(
            "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru)
    ]
    ai = interface.amber
    md_step = ai.build_md_step(
        length=0.1,
        constrain=test_constrains)
    test_params = AmberParameter(test_inpcrd, test_prmtop)
    test_job, test_md_egg = md_step.make_job(test_params)
    
    answer_pattern = r"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=MD_EnzyHTP
#SBATCH --partition=<fillthis>
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --account=<fillthis>
#SBATCH --export=NONE
(?:#SBATCH --exclude=.+)?

# Script generated by EnzyHTP [0-9]\.[0-9]\.[0-9] in [0-9]+-[0-9]+-[0-9]+ [0-9]+:[0-9]+:[0-9]+

source /home/shaoq1/bin/amber_env/amber22\.sh

pmemd\.cuda -O -i \./MD/amber_md_step_?[0-9]*\.in -o \./MD/amber_md_step\.out -p .*test/_interface/data//KE_07_R7_S\.prmtop -c .*test/_interface/data//KE_07_R7_S\.inpcrd -r \./MD/amber_md_step\.rst -ref .*test/_interface/data//KE_07_R7_S\.inpcrd -x \./MD/amber_md_step\.nc 
"""
    assert re.match(answer_pattern, test_job.sub_script_str)
    fs.safe_rmdir(md_step.work_dir)


def test_amber_md_step_try_merge_jobs(caplog):
    """test to make sure AmberMDStep.try_merge_jobs() works as expected."""
    with EnablePropagate(_LOGGER):
        ai = interface.amber
        md_step_1 = ai.build_md_step(length=0.1) # 300K, NPT by default
        md_step_2 = ai.build_md_step(length=0.1) # 300K, NPT by default
        md_step_3 = ai.build_md_step(length=0.1, core_type="cpu") # 300K, NPT by default
        test_inpcrd = f"{MM_DATA_DIR}/KE_07_R7_S.inpcrd"
        test_prmtop = f"{MM_DATA_DIR}/KE_07_R7_S.prmtop"
        test_params = AmberParameter(test_inpcrd, test_prmtop)
        test_job_1, test_md_egg_1 = md_step_1.make_job(test_params)
        test_job_2, test_md_egg_2 = md_step_2.make_job(test_md_egg_1)
        test_job_3, test_md_egg_3 = md_step_3.make_job(test_md_egg_2)
        merged_jobs = AmberMDStep.try_merge_jobs([test_job_1,test_job_2,test_job_3])
        assert len(merged_jobs) == 2
        for test, answer in zip(merged_jobs[0].mimo["commands"], [
                r'pmemd.cuda -O -i ./MD/amber_md_step_?[0-9]*.in -o ./MD/amber_md_step.out -p .*test/_interface/data//KE_07_R7_S.prmtop -c .*test/_interface/data//KE_07_R7_S.inpcrd -r ./MD/amber_md_step.rst -ref .*test/_interface/data//KE_07_R7_S.inpcrd -x ./MD/amber_md_step.nc ',
                r'pmemd.cuda -O -i ./MD/amber_md_step_?[0-9]*.in -o ./MD/amber_md_step.out -p .*test/_interface/data//KE_07_R7_S.prmtop -c ./MD/amber_md_step.rst -r ./MD/amber_md_step.rst -ref ./MD/amber_md_step.rst -x ./MD/amber_md_step.nc '
                ]):
            assert re.match(answer, test)
        assert len(merged_jobs[0].mimo["temp_mdin"]) == 2
        assert "Found md steps with same names!" in caplog.text
    fs.safe_rmdir(md_step_1.work_dir)


def test_amber_md_step_translate():
    """use a fake test result egg to test the function"""
    raise Exception("deprocated") # figure out a way to test it
    test_result_egg = AmberMDResultEgg(
        traj_path = "traj_path",
        traj_log_path = "traj_log_path",
        rst_path = "rst_path",
        prmtop_path = "prmtop_pat",
    )
    ai = interface.amber
    test_step = ai.build_md_step(length=0.1)
    result = test_step.translate(test_result_egg)

    assert result.traj_file == "traj_path"
    assert result.traj_log_file == "traj_log_path"
    assert result.last_frame_file == "rst_path"


@pytest.mark.accre
@pytest.mark.long
def test_amber_md_step_run():
    """test AmberMDStep().run().
    rely on check_md_error inside of run() for testing.
    NOTE need GPU interactive shell"""
    ai = interface.amber
    md_step = ai.build_md_step(length=0.01) # 300K, NPT by default
    test_inpcrd = f"{MM_DATA_DIR}/KE_07_R7_S.inpcrd"
    test_prmtop = f"{MM_DATA_DIR}/KE_07_R7_S.prmtop"
    test_params = AmberParameter(test_inpcrd, test_prmtop)
    test_md_result = md_step.run(test_params)
    
    assert test_md_result.last_frame_file == "./MD/amber_md_step.rst"
    fs.clean_temp_file_n_dir([
        test_md_result.last_frame_file,
        test_md_result.traj_log_file,
        "./MD",
    ])


def test_amber_md_step_check_md_error(caplog):
    """test this function using extract real example files from Amber runs"""
    test_traj = f"{MM_DATA_DIR}md_errors/equi_error_1.nc"
    test_mdout = f"{MM_DATA_DIR}md_errors/equi_error_1.out"
    test_prmtop = f"{MM_DATA_DIR}md_errors/equi_error_1.prmtop"
    test_clusterjob = ClusterJob(
        cluster=None,
        sub_script_str=None,
    )
    test_clusterjob.job_cluster_log = f"{MM_DATA_DIR}md_errors/equi_error_1.stdstream"

    ai = interface.amber
    md_step = ai.build_md_step(length=0.1) # 300K, NPT by default

    with EnablePropagate(_LOGGER):
        with pytest.raises(AmberMDError) as e:
            md_step.check_md_error(
                traj=test_traj,
                traj_log=test_mdout,
                prmtop=test_prmtop,
                stdstream_source=test_clusterjob,
            )
        assert "Amber MD didn't finish normally." in caplog.text
        assert "ERROR: Calculation halted.  Periodic box dimensions have changed too much from their initial values." in caplog.text


def test_get_restraintmask_bb_freeze():
    """test using example restraints"""
    ai = interface.amber
    test_pdb = f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb"
    test_stru = struct.PDBParser().get_structure(test_pdb)
    cons = create_backbone_freeze(test_stru)
    mask = ai.get_restraintmask(cons)
    assert mask == "'@C,CA,N'"


def test_get_restraintmask():
    """test using example restraints"""
    ai = interface.amber
    test_pdb = f"{MM_DATA_DIR}/index_tweak.pdb"
    test_stru = struct.PDBParser().get_structure(test_pdb)
    cons = create_cartesian_freeze(
        atoms=test_stru.atoms[1:20],
        topology=test_stru,)
    mask = ai.get_restraintmask(cons)
    assert mask == "'@2-20'"


def test_get_amber_mask(): # TODO
    """"""
    raise Exception("TODO")


def test_get_amber_atom_index(): # TODO
    """"""
    raise Exception("TODO")


def test_get_amber_index_mapper():
    """this dont work for the 1Q4T case.
    Update: 2024.11.6 seems it works for 1Q4T?"""
    ai = interface.amber
    test_pdb = f"{STRU_DATA_DIR}/1Q4T_ligand_test.pdb"
    test_stru = struct.PDBParser().get_structure(test_pdb)
    test_res_1 = test_stru.get("A.12")
    test_res_2 = test_stru.get("B.12")
    test_atom_1 = test_stru.get("D.371.N1A")

    index_mapper = ai.get_amber_index_mapper(test_stru)
    assert index_mapper["residue"][test_res_1][1] == 3
    assert index_mapper["residue"][test_res_2][1] == 143
    assert index_mapper["atom"][test_atom_1] == 4434


def test_parse_cons_to_raw_rs_dict():
    """test using KE and example cons"""
    test_pdb = f"{MM_DATA_DIR}/KE_07_R7_2_S.pdb"
    test_stru = struct.PDBParser().get_structure(test_pdb)
    test_cons = create_distance_constraint(
        "B.254.H2", "A.101.OE2", 2.4, test_stru)
    answer_raw_dict = {
        'ialtd': 0,
        'iat': [3976, 1579],
        'r1': 2.15,
        'r2': 2.35,
        'r3': 2.45,
        'r4': 2.65,
        'rk2': 200.0,
        'rk3': 200.0}
    ai = interface.amber
    test_raw_dict = ai._parse_cons_to_raw_rs_dict(test_cons)
    for k, v in test_raw_dict.items():
        if k in "r1 r2 r3 r4".split():
            assert np.isclose(v, answer_raw_dict[k], atol=1e-3)
        else:
            assert v == answer_raw_dict[k]


def test_reduce_path_in_mdin():
    """test using an path that is too long"""
    test_path = os.path.abspath("test_integration/equi_md_qm_spe_based_descriptor/work_dir/MD/rep_0/0.rs")
    ai = interface.amber

    result = ai._reduce_path_in_mdin(test_path)
    assert len(result) <= 70

def test_get_mmpbsa_energy():
    """test the function using old data generated by EnzyHTP 1.0"""
    ai = interface.amber
    stru_esm = ai.load_traj(
        prmtop_path=f"{MM_DATA_DIR}/mmpbsa_test_sol.prmtop",
        traj_path=f"{MM_DATA_DIR}/mmpbsa_test_sol_10f.nc",
        ref_pdb=f"{MM_DATA_DIR}/mmpbsa_test_sol.pdb"
    )
    ligand = select_stru(stru_esm.structure_0, "resi 290")
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {
            "partition" : "production",
            "account" : "yang_lab",}
    }
    answer = [
        1.6840000000016948, 
        -0.6804000000007093, 
        5.076199999997218, 
        0.3467999999985514, 
        -5.959400000001956, 
        -5.139099999993475, 
        -3.543699999996821, 
        -5.468999999998417, 
        -3.973000000000989, 
        4.527400000004005
    ]
    
    result = ai.get_mmpbgbsa_energy(
        stru_esm, ligand, work_dir=MM_WORK_DIR,
        cluster_job_config=cluster_job_config,
    )
    for r, a in zip(result, answer):
        assert np.isclose(r, a, atol=0.01)
    
    fs.clean_temp_file_n_dir([f"{MM_WORK_DIR}/mmpbsa_by_frames.csv"])

def test_run_ante_mmpbsa():
    """test the function"""
    ai = interface.amber
    dry_complex_out = f"{MM_WORK_DIR}/antemmpbsa_dc.prmtop"
    dry_receptor_out = f"{MM_WORK_DIR}/antemmpbsa_dr.prmtop"
    dry_ligand_out = f"{MM_WORK_DIR}/antemmpbsa_dl.prmtop"

    ai.run_ante_mmpbsa(
        complex_prmtop_in=f"{MM_DATA_DIR}/mmpbsa_test_sol.prmtop",
        dry_complex_out = dry_complex_out,
        dry_receptor_out = dry_receptor_out,
        dry_ligand_out = dry_ligand_out,
        strip_mask = ":WAT,Na+,Cl-",
        radii = "mbondi2",
        ligand_mask = ":290"
    )

    assert os.path.getsize(dry_complex_out)
    assert os.path.getsize(dry_receptor_out)
    assert os.path.getsize(dry_ligand_out)

    fs.clean_temp_file_n_dir([
        dry_complex_out,
        dry_receptor_out,
        dry_ligand_out,
    ])

def test_update_radii():
    """test the function"""
    ai = interface.amber
    test_out_path = f"{MM_WORK_DIR}/update_radii_test.prmtop"

    ai.update_radii(
        prmtop_path=f"{MM_DATA_DIR}/mmpbsa_test_sol.prmtop",
        out_path=test_out_path,
        radii="mbondi2"
    )

    assert os.path.getsize(test_out_path)

    fs.clean_temp_file_n_dir([test_out_path])

def test_make_mmpbgbsa_prmtop_files():
    """test the function"""
    ai = interface.amber
    stru_esm = ai.load_traj(
        prmtop_path=f"{MM_DATA_DIR}/mmpbsa_test_sol.prmtop",
        traj_path=f"{MM_DATA_DIR}/mmpbsa_test_sol_10f.nc",
        ref_pdb=f"{MM_DATA_DIR}/mmpbsa_test_sol.pdb"
    )
    test_dr_prmtop = f"{MM_WORK_DIR}/test_make_mmpbsa_prmtop_dr.prmtop"
    test_dl_prmtop = f"{MM_WORK_DIR}/test_make_mmpbsa_prmtop_dl.prmtop"
    test_dc_prmtop = f"{MM_WORK_DIR}/test_make_mmpbsa_prmtop_dc.prmtop"
    test_sc_prmtop = f"{MM_WORK_DIR}/test_make_mmpbsa_prmtop_sc.prmtop"

    ai.make_mmpbgbsa_prmtop_files(
        stru_esm=stru_esm,
        ligand_mask=":290",
        strip_mask=":WAT,Na+,Cl-",
        igb=5,
        temp_dr_prmtop = test_dr_prmtop,
        temp_dl_prmtop = test_dl_prmtop,
        temp_dc_prmtop = test_dc_prmtop,
        temp_sc_prmtop = test_sc_prmtop,
    )

    assert os.path.getsize(test_dr_prmtop)
    assert os.path.getsize(test_dl_prmtop)
    assert os.path.getsize(test_dc_prmtop)
    assert os.path.getsize(test_sc_prmtop)

    fs.clean_temp_file_n_dir([
      test_dr_prmtop,test_dl_prmtop,test_dc_prmtop,test_sc_prmtop,stru_esm.topology_source_file
    ])

def test_count_num_of_frames_traj():
    ai = interface.amber
    prmtop_path = f"{MM_DATA_DIR}/mmpbsa_test_sol.prmtop"
    traj = f"{MM_DATA_DIR}/mmpbsa_test_sol_10f.nc"
    result = ai.count_num_of_frames_traj(
        prmtop_path=prmtop_path,
        traj_path=traj,
    )
    assert result == 10
    
def test_run_mmpbsa():
    """test the function"""
    ai = interface.amber
    dry_complex = f"{MM_DATA_DIR}/mmpbsa_test_dc.prmtop"
    dry_receptor = f"{MM_DATA_DIR}/mmpbsa_test_dr.prmtop"
    dry_ligand = f"{MM_DATA_DIR}/mmpbsa_test_dl.prmtop"
    sol_complex = f"{MM_DATA_DIR}/mmpbsa_test_sol.prmtop"
    traj = f"{MM_DATA_DIR}/mmpbsa_test_sol_10f.nc"
    out_path = f"{MM_WORK_DIR}/mmpbsa_test.csv"
    cluster_job_config = ClusterJobConfig.from_dict({
        "cluster" : Accre(),
        "res_keywords" : {
        'core_type' : 'cpu',
        'nodes':'1',
        'node_cores' : '24',
        'job_name' : 'mmpbsa_EnzyHTP',
        'partition' : 'production',
        'mem_per_core' : '2G', # in GB.
        'walltime' : '24:00:00',
        'account' : 'yang_lab',
        }
    })

    ai.run_mmpbsa(
        dr_prmtop=dry_receptor,
        dl_prmtop=dry_ligand,
        dc_prmtop=dry_complex,
        sc_prmtop=sol_complex,
        traj_file=traj,
        out_path=out_path,
        solvent_model="pbsa",
        cluster_job_config=cluster_job_config,
        job_check_period = 5,
    )

    assert os.path.exists(out_path)
    fs.clean_temp_file_n_dir([out_path])

def test_run_mmpbsa_local():
    """test the function"""
    ai = interface.amber
    dry_complex = f"{MM_DATA_DIR}/mmpbsa_test_dc.prmtop"
    dry_receptor = f"{MM_DATA_DIR}/mmpbsa_test_dr.prmtop"
    dry_ligand = f"{MM_DATA_DIR}/mmpbsa_test_dl.prmtop"
    sol_complex = f"{MM_DATA_DIR}/mmpbsa_test_sol.prmtop"
    traj = f"{MM_DATA_DIR}/mmpbsa_test_sol_10f.nc"
    out_path = f"{MM_WORK_DIR}/mmpbsa_test.csv"

    ai.run_mmpbsa(
        dr_prmtop=dry_receptor,
        dl_prmtop=dry_ligand,
        dc_prmtop=dry_complex,
        sc_prmtop=sol_complex,
        traj_file=traj,
        out_path=out_path,
        solvent_model="pbsa",
        non_armer_cpu_num=10,
    )

    assert os.path.exists(out_path)
    fs.clean_temp_file_n_dir([out_path])

def test_parse_mmpbsa_result_not_by_frame():
    """test the function using an example data file"""
    ai = interface.amber
    test_dat_file = f"{MM_DATA_DIR}/mmpbsa_test.dat"

    result = ai.parse_mmpbsa_result(test_dat_file, by_frames=False)

    assert "gbsa" not in result
    assert result["pbsa"]["mean"]["DELTA TOTAL"] == -1.3130

def test_parse_mmpbsa_result_by_frames():
    """test the function using an example data file. Using manually confirmed data as reference."""
    ai = interface.amber
    test_dat_file = f"{MM_DATA_DIR}/mmpbsa_test.csv"

    result = ai.parse_mmpbsa_result(test_dat_file, by_frames=True)

    assert np.isclose(result["pbsa"]["DELTA TOTAL"].mean(), -1.3130, atol=0.01)
    assert np.isclose(result["gbsa"]["DELTA TOTAL"].mean(), -12.234, atol=0.01)
