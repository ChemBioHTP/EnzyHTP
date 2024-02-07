#! python3
# -*- encoding: utf-8 -*-
'''
Test the enzy_htp.workflow.workflow module.

@File    :   test_workflow.py
@Created :   2024/02/07 01:30
@Author  :   Zhong, Yinjie
@Version :   1.0
@Contact :   yinjie.zhong@vanderbilt.edu
'''

# Here put the import lib.
import os
import pytest

from enzy_htp import config
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core.general import EnablePropagate

from enzy_htp.structure import structure
from enzy_htp.workflow.workflow import WorkFlow, ProcedureUnit

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data/"
WORK_DIR = f"{CURR_DIR}/work_dir/"
config["system.SCRATCH_DIR"] = WORK_DIR

def test_procedure_unit_self_inspection_read_pdb_7si9(caplog):
    '''A test for reading pdb file.'''
    pdb_path = f'{DATA_DIR}7si9_rm_water_disconnected.pdb'
    unit_dict = {
        "api" : "read_pdb",
        "return" : "read_pdb_0",
        "args" : {
            "path" : pdb_path,
        }
    }
    procedure_unit = ProcedureUnit.from_dict(unit_dict=unit_dict)
    stru: structure = procedure_unit.execute()
    assert stru.num_residues
    
def test_procedure_unit_self_inspection_unrecognized_api(caplog):
    '''A test for unrecognized API.'''
    pdb_path = f'{DATA_DIR}7si9_rm_water_disconnected.pdb'
    unit_dict = {
        "api" : "read",
        "return" : "read_something",
        "args" : {
            "path" : pdb_path,
        }
    }
    with pytest.raises(KeyError) as e:    # The value error is expected to be raised when unexpected Hydrogen atom(s) is detected.
            with EnablePropagate(_LOGGER):
                procedure_unit = ProcedureUnit.from_dict(unit_dict=unit_dict)
    assert 'any api' in caplog.text.lower()

def test_procedure_unit_self_inspection_missing_arg(caplog):
    '''A test for missing required arguments.'''
    unit_dict = {
        "api" : "read_pdb",
        "return" : "read_pdb_0",
        "args" : {
            # "path" : "xxx",
        }
    }
    with pytest.raises(ValueError) as e:    # The value error is expected to be raised when unexpected Hydrogen atom(s) is detected.
            with EnablePropagate(_LOGGER):
                procedure_unit = ProcedureUnit.from_dict(unit_dict=unit_dict)
    assert 'miss' in caplog.text.lower()
    
def test_procedure_unit_self_inspection_unexpected_type(caplog):
    '''A test for unexpected argument types.'''
    unit_dict = {
        "api" : "read_pdb",
        "return" : "read_pdb_0",
        "args" : {
            "path" : 1,
            "model": 'xxx',
        }
    }
    with pytest.raises(ValueError) as e:    # The value error is expected to be raised when unexpected Hydrogen atom(s) is detected.
            with EnablePropagate(_LOGGER):
                procedure_unit = ProcedureUnit.from_dict(unit_dict=unit_dict)
    assert 'expect' in caplog.text.lower()