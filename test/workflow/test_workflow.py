#! python3
# -*- encoding: utf-8 -*-
'''
Test the enzy_htp.task_assembler.workflow module.

@File    :   test_workflow.py
@Created :   2024/02/07 01:30
@Author  :   Zhong, Yinjie
@Version :   1.0
@Contact :   yinjie.zhong@vanderbilt.edu
'''

# Here put the import lib.
import os
import pytest
import logging

from typing import Any
from datetime import datetime
from random import randint

from enzy_htp import config
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core.clusters.accre import Accre
from enzy_htp.core.general import EnablePropagate

from enzy_htp.structure import Structure
from enzy_htp.chemical.level_of_theory import QMLevelOfTheory

from enzy_htp.workflow import ExecutionEntity, WorkFlow, WorkUnit, GeneralWorkUnit
from enzy_htp.workflow import SCIENCE_API_MAPPER, StatusCode

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/data"
WORK_DIR = f"{CURR_DIR}/work_dir"
config["system.SCRATCH_DIR"] = WORK_DIR

class PseudoAPI():
    """This class acts as a container to save some static pseudo SCIENCE APIs for test."""

    @staticmethod
    def test_kwargs(x: str, **kwargs):
        '''This function is to test kwargs inspection only.'''
        result = x
        for key, value in kwargs.items():
            result += f'\n {key} = {value}'
        return result
    
    @staticmethod
    def log_atom_number_difference(original_stru: Structure, mutant_stru: Structure):
        """This function is to print the difference between WildType and Mutant in atom numbers."""
        _LOGGER.info(f"The original structure contains {original_stru.num_atoms} atoms.")
        _LOGGER.info(f"The mutant structure contains {mutant_stru.num_atoms} atoms.")
        return {"original": original_stru.num_atoms, "mutant": mutant_stru.num_atoms}
    
    @staticmethod
    def union_values(**kwargs) -> str:
        """This function prints out all the keyword argument values and return the printed value."""
        unioned_value = " | ".join(map(str, kwargs.values()))
        _LOGGER.info(f"Unioned Value is: {unioned_value}")
        return unioned_value
    
    @staticmethod
    def print_value(value) -> None:
        print(value)
        return
    
    @staticmethod
    def checkpoint_fruit_producer(fruit_name: str) -> None:
        """This function produce input value for checkpoint functions.
        
        Args:
            fruit_name (str): The name of the fruit you want to produce.
        
        Returns:
            A generated string value containing the `fruit_name` you enter.
        """
        brands = ["Mutexa", "Enzeva", "Minerva", "Molgraph"]
        current_time = datetime.now()
        result = f"{brands[randint(0, len(brands)-1)]}'s {fruit_name}, produced at {current_time.strftime('%Y-%m-%d_%H:%M')}."
        return result
        
    @staticmethod
    def checkpoint_apple_eater_error(test_str: str) -> None:
        """This function act as a checkpoint to test if the reload (continue computing) functions normally.
        
        Args:
            test_str (str): A string value as input for checkpoint.

        Raises:
            ValueError: If the `test_str` does not contain `apple`.
        """
        _LOGGER.info(f"I received {test_str}...")
        if ("apple" in test_str):
            return
        else:
            raise ValueError("I need apples!")

PSEUDO_API_MAPPER = {
    "test_kwargs": PseudoAPI.test_kwargs,
    "log_difference": PseudoAPI.log_atom_number_difference,
    "print": PseudoAPI.print_value,
    "union_values": PseudoAPI.union_values,
    "checkpoint_producer": PseudoAPI.checkpoint_fruit_producer,
    "checkpoint_error": PseudoAPI.checkpoint_apple_eater_error
}

SCIENCE_API_MAPPER.update(PSEUDO_API_MAPPER)

def test_workunit_read_pdb_7si9(caplog):
    '''A test for reading pdb file.'''
    pdb_path = f'{DATA_DIR}/7si9_rm_water_disconnected.pdb'
    unit_dict = {
        "api" : "read_pdb",
        "store_as" : "read_pdb_0",
        "args" : {
            "path" : pdb_path,
        }
    }
    workunit = WorkUnit.from_dict(unit_dict=unit_dict, debug=True)
    key, stru = workunit.execute()
    assert stru.num_residues
    
def test_workunit_self_inspection_unrecognized_api(caplog):
    '''A test for unrecognized API.'''
    pdb_path = f'{DATA_DIR}/7si9_rm_water_disconnected.pdb'
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

def test_workunit_self_inspection_kwargs_pseudo_api(caplog):
    '''Check if kwargs can be successfully inspected.'''
    unit_dict = {
        "api" : "test_kwargs",
        "store_as" : "assign_mutant_0",
        "args" : {
            "x" : "Let's Rock!",
            "pattern" : "H41Y, M165S",
            "chain_sync_list" : ["A", "B"],
            "chain_index_mapper" : {"A" : 0, "B" : 100},
        }
    }
    workunit = WorkUnit.from_dict(unit_dict=unit_dict, debug=True)
    key, result = workunit.execute()
    print(result)
    assert 'chain' in result

def test_workflow_store_as_label(caplog):
    '''Test whether the updated `intermediate_data_mapper` labeling strategy 
    (i.e. use the same label for instances with same memory address) works.'''
    json_filepath = f'{DATA_DIR}/workflow_7si9_store_as_label.json'
    workflow = WorkFlow.from_json_filepath(json_filepath=json_filepath)
    _LOGGER.level = logging.DEBUG
    with EnablePropagate(_LOGGER):
        workflow.execute()
    fs.safe_rmdir(WORK_DIR)
    print(workflow.intermediate_data_mapper)
    assert workflow.intermediate_data_mapper['structure']

def test_workflow_initialization_from_json_filepath(caplog):
    '''Test initializing workflow from a simple json file.
    This test may have errors reported, but they are expected.
    '''
    json_filepath = f'{DATA_DIR}/workflow_7si9_initialization.json'
    workflow = WorkFlow.from_json_filepath(json_filepath=json_filepath)
    with EnablePropagate(_LOGGER):
        workflow.execute()
    fs.safe_rmdir(WORK_DIR)
    # preparation.protonate_stru does not have return value.
    # print(workflow.intermediate_data_mapper)
    assert 'success' in caplog.text
    return

def test_general_loopworkunit_pseudo_api(caplog):
    """Test initializing and executing workflow with Loop(s) to test if LoopWorkUnit works properly."""
    json_filepath = f'{DATA_DIR}/general_7si9_loopworkunit_pseudo_api.json'
    # _LOGGER.level = logging.DEBUG
    with EnablePropagate(_LOGGER):
        general = GeneralWorkUnit.from_json_filepath(json_filepath=json_filepath, working_directory=WORK_DIR, overwrite_database=True, debug=True)
        return_key, return_value = general.execute()
    fs.safe_rmdir(WORK_DIR)
    assert 'original structure' in caplog.text
    assert 'mutant structure' in caplog.text
    return

def test_general_layer_variable(caplog):
    """Test initializing and executing general with all the variables 
    with all user-defined arguments are set in the general layer."""
    json_filepath = f'{DATA_DIR}/general_7si9_general_layer_variable.json'
    # _LOGGER.level = logging.DEBUG
    with EnablePropagate(_LOGGER):
        general = GeneralWorkUnit.from_json_filepath(json_filepath=json_filepath, working_directory=WORK_DIR, overwrite_database=True, debug=True)
        return_key, return_value = general.execute()
    fs.safe_rmdir(WORK_DIR)
    assert 'original structure' in caplog.text
    assert 'mutant structure' in caplog.text
    return

def test_general_multi_loop_datum(caplog):
    """Test initializing and executing general with multi layer loops 
    to check if loop_datum can be passed correctly between multi layers."""
    json_filepath = f'{DATA_DIR}/general_multi_loop_datum.json'
    # _LOGGER.level = logging.DEBUG
    with EnablePropagate(_LOGGER):
        general = GeneralWorkUnit.from_json_filepath(json_filepath=json_filepath, working_directory=WORK_DIR, overwrite_database=True, debug=True)
        return_key, return_value = general.execute()
    fs.safe_rmdir(WORK_DIR)
    assert 'Lisa | Green Pepper' in caplog.text
    return

def test_general_locate_workunit(caplog):
    """Test locating a certain workunit after the execution of the general."""
    json_filepath = f'{DATA_DIR}/general_7si9_general_layer_variable.json'
    general = GeneralWorkUnit.from_json_filepath(json_filepath=json_filepath, working_directory=WORK_DIR, overwrite_database=True, debug=True)
    return_key, return_value = general.execute()
    locator = GeneralWorkUnit.get_locator("mutate_stru@4:loop_2:0")
    target = general.locate(locator)
    fs.safe_rmdir(WORK_DIR)
    target_info = {
        "api_key": target.api_key,
        "locator": target.locator,
        "status": target.status
    }
    assert target_info["api_key"] == "mutate_stru"
    assert target_info["status"] == 0
    return

# TODO (Zhong): A test for database.

def test_general_read_architecture(caplog):
    """Test reading the schema of the configuration json/unit_dict."""
    json_filepath = f'{DATA_DIR}/general_7si9_general_layer_variable.json'
    general = GeneralWorkUnit.from_json_filepath(json_filepath=json_filepath, working_directory=WORK_DIR, overwrite_database=True, debug=True)
    # print(general.architecture)
    assert "        - mutate_stru" in general.architecture
    fs.safe_rmdir(WORK_DIR)
    return

def test_general_dump_and_load(caplog):
    """Test dumping pickle when facing error, and then loading pickle to continue computing."""
    json_filepath_error = f"{DATA_DIR}/general_7si9_reload_checkpoint_error.json"
    general = GeneralWorkUnit.from_json_filepath(json_filepath=json_filepath_error, working_directory=WORK_DIR, save_snapshot=True, overwrite_database=True)
    general.execute()

    # print("Finished first execution.")
    # print(general.architecture)
    # input()

    pickle_filepath = general.latest_pickle_filepath

    loop_body = general.locate_by_identifier(identifier="workflow@5:loop_1")
    mutate_stru = general.locate_by_identifier(identifier="mutate_stru@5:loop_1:0")
    checkpoint = general.locate_by_identifier(identifier="checkpoint_error@5:loop_2:1")
    assert loop_body.status == StatusCode.EXIT_WITH_ERROR_IN_INNER_UNITS
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.EXIT_WITH_ERROR

    json_filepath_pass = f"{DATA_DIR}/general_7si9_reload_checkpoint_pass.json"

    reloaded_general = GeneralWorkUnit.load_snapshot_file(filepath=pickle_filepath)
    reloaded_general.reload_json_filepath(json_filepath=json_filepath_pass)
    # print("Finished reloading.")
    # input()

    loop_body = reloaded_general.locate(ExecutionEntity.get_locator(identifier="workflow@5:loop_1"))
    mutate_stru = reloaded_general.locate(ExecutionEntity.get_locator(identifier="mutate_stru@5:loop_1:0"))
    checkpoint = reloaded_general.locate(ExecutionEntity.get_locator(identifier="checkpoint_error@5:loop_2:1"))
    assert loop_body.status == StatusCode.SUSPECIOUS_UPDATES
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.SUSPECIOUS_UPDATES

    reloaded_general.execute()
    # print("Finished second execution.")
    # input()

    assert loop_body.status == StatusCode.EXIT_OK
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.EXIT_OK

    fs.safe_rmdir(WORK_DIR)
    return

def test_general_reload_pass(caplog):
    """Test reloading and continue computing when facing error."""
    json_filepath_error = f"{DATA_DIR}/general_7si9_reload_checkpoint_error.json"
    general = GeneralWorkUnit.from_json_filepath(json_filepath=json_filepath_error, working_directory=WORK_DIR, overwrite_database=True)
    general.execute()

    loop_body = general.locate(ExecutionEntity.get_locator(identifier="workflow@5:loop_1"))
    mutate_stru = general.locate(ExecutionEntity.get_locator(identifier="mutate_stru@5:loop_1:0"))
    checkpoint = general.locate(ExecutionEntity.get_locator(identifier="checkpoint_error@5:loop_2:1"))

    assert loop_body.status == StatusCode.EXIT_WITH_ERROR_IN_INNER_UNITS
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.EXIT_WITH_ERROR

    json_filepath_pass = f"{DATA_DIR}/general_7si9_reload_checkpoint_pass.json"
    general.reload_json_filepath(json_filepath=json_filepath_pass)
    assert loop_body.status == StatusCode.SUSPECIOUS_UPDATES
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.SUSPECIOUS_UPDATES

    general.execute()
    assert loop_body.status == StatusCode.EXIT_OK
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.EXIT_OK

    fs.safe_rmdir(WORK_DIR)
    return

def test_general_reload_change_varname(caplog):
    """Test reloading and test if changed varname can be detecetd and raise error."""
    json_filepath_error = f"{DATA_DIR}/general_7si9_reload_checkpoint_error.json"
    general = GeneralWorkUnit.from_json_filepath(json_filepath=json_filepath_error, working_directory=WORK_DIR, overwrite_database=True)
    general.execute()

    loop_body = general.locate(ExecutionEntity.get_locator(identifier="workflow@5:loop_1"))
    mutate_stru = general.locate(ExecutionEntity.get_locator(identifier="mutate_stru@5:loop_1:0"))
    checkpoint = general.locate(ExecutionEntity.get_locator(identifier="checkpoint_error@5:loop_2:1"))

    assert loop_body.status == StatusCode.EXIT_WITH_ERROR_IN_INNER_UNITS
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.EXIT_WITH_ERROR

    json_filepath_pass = f"{DATA_DIR}/general_7si9_reload_checkpoint_change_varname.json"
    general.reload_json_filepath(json_filepath=json_filepath_pass)
    assert loop_body.status == StatusCode.SUSPECIOUS_UPDATES
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.SUSPECIOUS_UPDATES

    general.execute()
    assert loop_body.status == StatusCode.EXIT_OK
    assert mutate_stru.status == StatusCode.EXIT_OK
    assert checkpoint.status == StatusCode.EXIT_OK

    fs.safe_rmdir(WORK_DIR)
    return

def test_general_execute_full_workflow(caplog):
    """Execute the full workflow."""
    cluster = Accre()
    md_cluster_job_config = {
        "cluster" : cluster,
        "res_keywords" : {
            "account" : "csb_gpu_acc",
            "partition" : "a6000x4"
        }
    }
    qm_cluster_job_config = {
        "cluster" : cluster,
        "res_keywords" : {
            "account" : "yang_lab_csb",
            "partition" : "production",
            'walltime' : '1-00:00:00',
        }
    }
    qm_level_of_theory = QMLevelOfTheory(
        basis_set="3-21G",
        method="HF",
        solvent="water",
        solv_method="SMD",
    )
    data_mapper_for_init = {
        "cluster": cluster,
        "md_cluster_job_config": md_cluster_job_config,
        "qm_cluster_job_config": qm_cluster_job_config,
        "qm_level_of_theory": qm_level_of_theory,
    }

    json_filepath = f"{DATA_DIR}/general_7si9_general_full.json"
    general = GeneralWorkUnit.from_json_filepath(
        json_filepath=json_filepath, working_directory=WORK_DIR, 
        overwrite_database=True, data_mapper_for_init=data_mapper_for_init)
    pass
    