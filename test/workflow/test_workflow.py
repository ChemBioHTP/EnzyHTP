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
from enzy_htp.workflow import WorkFlow, WorkUnit
from enzy_htp.workflow.config import SCIENCE_API_MAPPER

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data/"
WORK_DIR = f"{CURR_DIR}/work_dir/"
config["system.SCRATCH_DIR"] = WORK_DIR


def test_workflow_initialization_from_json_filepath_simple(caplog):
    '''Test initializing workflow from a simple json file.'''
    json_filepath = f'{DATA_DIR}workflow_7si9_init_simple.json'
    workflow = WorkFlow.from_json_filepath(json_filepath=json_filepath)
    with EnablePropagate(_LOGGER):
        workflow.execute()
    fs.safe_rmdir(WORK_DIR)
    # preparation.protonate_stru does not have return value.
    assert workflow.intermediate_data_mapper['protonate_stru_0']

def test_workflow_initialization_from_json_filepath_typing(caplog):
    json_filepath = f'{DATA_DIR}workflow_7si9_init_typing.json'
    workflow = WorkFlow.from_json_filepath(json_filepath=json_filepath)
    with EnablePropagate(_LOGGER):
        workflow.execute()
    fs.safe_rmdir(WORK_DIR)
    # [02:46 EnzyHTP <listcomp>] ERROR Receiving argument `chain_sync_list` in unexpected type <class 'list'> while typing.List[tuple] is expected (when initializing `assign_mutant`).
    # [02:46 EnzyHTP <listcomp>] ERROR Receiving argument `chain_index_mapper` in unexpected type <class 'dict'> while typing.Dict[str, int] is expected (when initializing `assign_mutant`).
    # [02:46 EnzyHTP execute] ERROR The initialized WorkFlow has not yet passed the self-inspection, so it is not allowed to be executed.
    assert 'success' in caplog.text

def test_workunit_read_pdb_7si9(caplog):
    '''A test for reading pdb file.'''
    pdb_path = f'{DATA_DIR}7si9_rm_water_disconnected.pdb'
    unit_dict = {
        "api" : "read_pdb",
        "store_as" : "read_pdb_0",
        "args" : {
            "path" : pdb_path,
        }
    }
    workunit = WorkUnit.from_dict(unit_dict=unit_dict, debug=True)
    key, stru = workunit.execute()
    print(stru.num_residues)
    assert stru.num_residues

def test_workunit_self_inspection_kwargs(caplog):
    '''Check if kwargs can be successfully inspected.'''
    from inspect import signature

    unit_dict = {
        "api" : "test_kwargs",
        "store_as" : "assign_mutant_0",
        "args" : {
            "x" : "Let's Rock!",
            "pattern" : "xxx",
            "chain_sync_list" : ["A", "B"],
            "chain_index_mapper" : {"A" : 0, "B" : 100},
        }
    }
    workunit = WorkUnit.from_dict(unit_dict=unit_dict, debug=True)
    key, result = workunit.execute()
    assert 'chain' in result
    
def test_workunit_self_inspection_unrecognized_api(caplog):
    '''A test for unrecognized API.'''
    pdb_path = f'{DATA_DIR}7si9_rm_water_disconnected.pdb'
    unit_dict = {
        "api" : "read",
        "store_as" : "read_something",
        "args" : {
            "path" : pdb_path,
        }
    }
    with pytest.raises(KeyError) as e:    # The value error is expected to be raised when unexpected Hydrogen atom(s) is detected.
            with EnablePropagate(_LOGGER):
                workunit = WorkUnit.from_dict(unit_dict=unit_dict, debug=True)
    assert 'map' in caplog.text.lower()

def test_workunit_self_inspection_missing_arg(caplog):
    '''A test for missing required arguments.'''
    unit_dict = {
        "api" : "read_pdb",
        "store_as" : "read_pdb_0",
        "args" : {
            # "path" : "xxx",
        }
    }
    with pytest.raises(ValueError) as e:    # The value error is expected to be raised when unexpected Hydrogen atom(s) is detected.
            with EnablePropagate(_LOGGER):
                workunit = WorkUnit.from_dict(unit_dict=unit_dict, debug=True)
    assert 'miss' in caplog.text.lower()
    
def test_workunit_self_inspection_unexpected_type(caplog):
    '''A test for unexpected argument types.'''
    unit_dict = {
        "api" : "read_pdb",
        "store_as" : "read_pdb_0",
        "args" : {
            "path" : 1,
            "model": 'xxx',
        }
    }
    with pytest.raises(ValueError) as e:    # The value error is expected to be raised when unexpected Hydrogen atom(s) is detected.
            with EnablePropagate(_LOGGER):
                workunit = WorkUnit.from_dict(unit_dict=unit_dict, debug=True)
    assert 'expect' in caplog.text.lower()