#! python3
# -*- encoding: utf-8 -*-
"""
Workflow Interpreter.

TODO (Zhong): Typing type annotation inspection. (Completed 2024-02-16 23:16 UTC-6)
TODO (Zhong): Annotation inspection in intermediate_data_mapper. (Completed 2024-03-01 01:00 UTC-6)
TODO (Zhong): APIs entitled `loop`. (Completed 2024-02-26 20:01 UTC-6)
TODO (Zhong): How to handle intermediate data transfer between parent and child workflows. Child workflow should be able to read the data in parent flow. (Completed 2024-02-28 06:06 UTC-6)
TODO (Zhong): Continue computing. (Completed 2024-03-20 01:15 UTC-5)
TODO (Zhong): How to indicate error from inner units. (Completed 2024-03-13 19:00 UTC-5)
TODO (Zhong): Cluster Batch.
TODO (Zhong): Data Output to excel (enzy_htp.core.file_system.write_data_to_excel, Completed 2024-05-12 19:50 UTC-5). 

@File    :   workflow.py
@Created :   2024/02/05 22:53
@Author  :   Zhong, Yinjie
@Contact :   yinjie.zhong@vanderbilt.edu
"""

# Here put the import lib.
from __future__ import annotations
from os import getcwd, chdir, path
from io import TextIOWrapper, FileIO
from time import sleep
from inspect import _empty, Parameter, signature
from json import load, loads
from typing import Any, Dict, List, Type, Union, Tuple, get_args
from collections.abc import Iterable, Callable
from datetime import datetime
from sqlalchemy.orm.session import Session
import pickle

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.file_system import safe_mkdir

from .config import (
    SCIENCE_API_MAPPER, 
    JsonConfigVariableFormatter as VarFormatter,
    StatusCode,
    Placeholder,
)

from .database import (
    SqliteSqlalchemy,
    ExecutionEntity as EntityORM
)

BASE_DIRECTORY = getcwd()   # Base working directory. 
                            # (Not the current directory where this file is located, but the directory where the program is triggered.)
DEFAULT_SQLITE_FILENAME = "database.db" # Default database filename.

# In the Json file used to define the workflow, 
# keys in the mapper are defined by the developer,
# which are not user-definable.
WORKING_DIRECTORY_KEY = "working_directory"     # Key indicating the working directory in the data mapper.
WORKUNIT_API_NAME_KEY = "api"                   # Key indicating the API name.
WORKUNIT_RETURN_VALUE_KEY = "store_as"          # Key indicating the key where the return value of the workunit is stored.
WORKUNIT_ARGUMENT_LIST_KEY = "args"             # Key indicates the arguments to be passed into the API.
CONTROL_BODY_WORKUNITS_LABEL = "workunits"      # Key indicating the list of workunits contained in the body in a unit playing a control role.
WORKFLOW_API_KEY = "workflow"                   # The value of `api_key` attribute of workflow instance.
LOOP_API_KEY = "loop"                           # The key of a loop workunit.
LOOP_ITERABLE_DATA_LABEL = "loop_data"          # Key indicating the object to iterate over in a LoopWorkUnit.
LOOP_BODY_DATUM_LABEL = "loop_datum_varname"    # Key indicating the element iterated from ITERABLE_DATA in each workunit of the loop body.
CLUSTER_BATCH_API_KEY = "cluster_batch"         # The key of a cluster batch workunit.
BATCH_ITERABLE_DATA_LABEL = "batch_data"        # Key indicating the object to iterate over in a ClusterBatchWorkUnit.
BATCH_BODY_DATUM_LABEL = "batch_datum_varname"  # Key indicating the element iterated from ITERABLE_DATA in each workunit of the batch body.
DEFAULT_CLUSTER_JOB_CAPABILITY = 5              # How much jobs can be submitted to the cluster before reaching its capability?
GENERAL_API_KEY = "general"                     # The key of a general workunit.

# API Keys for flow control, such as `loop`, `general`, etc.
CONTROL_API_KEYS = [
    LOOP_API_KEY,
    CLUSTER_BATCH_API_KEY,
    GENERAL_API_KEY,
]

# Some string values used globally are defined here.
FAILED_INITIALIZATION_ERROR_MSG = "The initialized work has not yet passed the self-inspection, so it is not allowed to be executed!"
SUCCESS_INITIALIZATION_MSG = "The initialized work has successfully passed the self-inspection. Now proceeding to execution."

#region Base entities.

class ExecutionEntity:
    """
    A base class representing an entity in an execution context, such as a workflow or workunit.

    This class provides a shared foundation for both WorkUnit and WorkFlow classes, encapsulating
    common attributes and functionalities like identification and executability check. It also
    facilitates the recursive search of entities based on their unique identifiers, supporting 
    complex execution structures in workflows.

    Attributes:
        api_key: str    # A string value indicating the role/function of the current ExecutionEntity instance.
        status (int): A integer indicating the execution status of the current ExecutionEntity instance.
                    The value of `status` is according to the values in `StatusCode` class.
        locator (List[int]): A list of the hierarchical indexes to locate the current ExecutionEntity instance.
                            If the current instance is at the outermost level, this value is an empty list.
        debug (bool): Flag indicating whether debugging mode is enabled.
                    If True, the error is raised everytime it is catched.
        child_execution_entities (list): Child execution entities of current instance.
                                        e.g. `sub_workflows` in `LoopWorkUnit`, `workunits` in `WorkFlow`, etc.
        working_directory (str): Indicates where to run the job.
        sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
        identifier: A unique identifier for the instance, used to track and reference the entity.
        is_executable: Indicates whether the instance has passed self-inspection and is ready for execution.

    Methods:
        locate: Recursively searches for an ExecutionEntity instance within the workflow or process
                structure based on the locator. Useful for locating specific components within
                a potentially complex workflow.
    """
    api_key: str    # A string value indicating the role/function of the current ExecutionEntity instance.
    _status: int    # A integer indicating the execution status of the current ExecutionEntity instance.
    locator: list   # A list of the hierarchical indexes to locate the current ExecutionEntity instance.
    debug: bool     # Indicates whether to run in debug mode.
    child_execution_entities: list  # Child execution entities of current instance.
    database_session: Session   # The database connection session.
    working_directory: str      # Indicates where to run the job.
    sqlite_filepath: str        # The filepath to your sqlite database file for persistent storage. 

    def __init__(self, working_directory: str = BASE_DIRECTORY, debug: bool = False, locator: List[int] = list(), sqlite_filepath: str = str()):
        """Initialize an ExecutionEntity instance.
        
        Args:
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            debug (bool, optional): Indicates whether to run in debug mode.
            locator (List[int], optional): A list of the hierarchical indexes to locate the current workflow.
                                If the current workflow is at the outermost level, this value is an empty list.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
        """
        self.working_directory = working_directory
        self.debug = debug
        self.locator = locator
        self.sqlite_filepath = sqlite_filepath
        self._status = StatusCode.CREATED
        self.api_key = str()
        self.child_execution_entities = None
        if (sqlite_filepath):
            self.sqlite_filepath = sqlite_filepath

    def create_database_session(self, overwrite_database: bool = False) -> Session:
        """Create database session with given working directory and sqlite_filename.
        
        Args:
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            overwrite_database: If overwrite the database file when it exists. Default False.

        Returns:
            The database connection session to save entity status as persistent storage.
        """
        session = SqliteSqlalchemy(self.sqlite_filepath, overwrite_database).session
        return session
    
    @property
    def status(self) -> int:
        return self._status
    
    @status.setter
    def status(self, value: int):
        self._status = value
        if (value not in StatusCode.unexecuted_statuses):
            self.synchronize_execution_status()
    
    def synchronize_execution_status(self) -> None:
        """Synchronize the status of the current entity to database. Called everytime before executing the WorkUnit or WorkFlow instances.
        
        Note: Don't call it within initialization!

        TODO (Zhong): 也许一个根据 Walltime 时间进行重生的逻辑可以放在这个卡口。Maybe a logic for respawning based on Walltime time could go in this place.
        """
        if (not self.is_executable):
            return
        if (not self.sqlite_filepath):
            _LOGGER.warning(f"No SQLite filepath specified, so the execution status of {self.identifier} is not synchronized!")
            _LOGGER.warning("If this warning appears during initialization, ignore it.")
            return
        database_session = self.create_database_session()
        if (database_session):
            entity = EntityORM(self.identifier, self.status)
            entity.insert_or_update(database_session)
            database_session.close()
        else:
            _LOGGER.warning(f"Database session doesn't exist in {self.identifier}!")
            return
    
    @property
    def is_executable(self) -> bool:
        """Flag indicating whether the current instance passed self-inspection."""
        if self.status in StatusCode.unexecutable_statuses:
            return False
        else:
            return True

    @property
    def identifier(self):
        """The unique identifier of current instance, collapsed from `api_key` and `locator`."""
        identifier = f"{self.api_key}@{':'.join(map(str, self.locator))}"
        return identifier
    
    @staticmethod
    def get_locator(identifier: str) -> list:
        """Parse the identifier string of an ExecutionEntity instance and return the locator list.
        
        Args:
            identifier (str): A unique identifier for the instance, used to track and reference the entity.

        Returns:
            The list of locators corresponding to the identifier.
        """
        locator_section = identifier.split("@")[-1]
        raw_locator = locator_section.split(":")
        locator = []
        for item in raw_locator:
            try:
                # Try parsing item to int.
                locator.append(int(item))
            except ValueError:
                # Keep original item if it's not an int.
                locator.append(item)
        return locator

    def locate(self, relative_locator: list) -> ExecutionEntity:
        """Recursive locating of ExecutionEntity instance based on their locators.

        Note 1: This method will return None for an unexecuted ExecutionEntity instance that is inside a ControlWorkUnit instance other than GeneralWorkUnit.
        Note 2: This method is overwritten in ControlWorkUnit and GeneralWorkUnit classes due to different behaviors.
        
        Args:
            relative_locator (list): A locator that locates the target instance starting from the current instance.
                                    For example, if the locator of target instance is [4, loop_0, 1, 2] and the locator of current instance is [4, loop_0],
                                    then the value of `relative_locator` should be [1, 2].

        Returns:
            The located ExecutionEntity instance. None if nothing is located.
        """
        if len(relative_locator) == 0:
            return self
        else:
            try:
                located_instance: ExecutionEntity = self.child_execution_entities[relative_locator[0]]
                return located_instance.locate(relative_locator[1:])
            except:
                return None
    
    def locate_by_identifier(self, identifier: str) -> ExecutionEntity:
        """Locate an ExecutionEntity instance by identifier.
        
        Args:
            identifier (str): A unique identifier for the instance, used to track and reference the entity.

        Returns:
            The located ExecutionEntity instance. None if nothing is located.
        """
        locator = __class__.get_locator(identifier=identifier)
        entity = self.locate(locator)
        if (entity):
            if (entity.identifier != identifier):
                return None
        return entity

    def dump_snapshot_file(self, save_directory: str = str()) -> str:
        """
        Saves the current state of the GeneralWorkUnit instance to a snapshot file.

        This method serializes the entire state of the GeneralWorkUnit instance and
        saves it to a file in the specified directory. The filename is generated
        based on the current timestamp.

        Note:
            Once you dump your snapshot file, the database session will be dumped.

            If the `save_directory` value entered is not a directory, 
            or if saving a file to this directory fails, the save directory would be 
            switched to the `working_directory` and a Warning is prompted; 
            if the file cannot be saved to the working directory either, an exception would be raised.

        Args:
            save_directory (str, optional): The directory where the snapshot file will be saved.
                                        Defaults to `working_directory`.

        Returns:
            str: The path to the created snapshot file.
        """
        dir_to_save = self.working_directory
        if (path.exists(save_directory)):
            dir_to_save = save_directory
        elif (save_directory):
            _LOGGER.warning(f"The `save_directory` value '{save_directory}' you specified is not a directory, so change the save address to working_directory.")
        else:
            _LOGGER.info("No specified directory to save the Pickle file. Save it to the working directory.")
        
        current_time = datetime.now()
        filename = f"snapshot_{self.api_key}_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}.pickle"
        filepath = path.join(dir_to_save, filename)
        try:
            with open(file=filepath, mode="wb") as fobj:
                pickle.dump(self, file=fobj)
        except Exception as exc:
            if dir_to_save != self.working_directory:
                _LOGGER.warning(f"Unable to save file to '{dir_to_save}'. Using your working directory '{self.working_directory}' instead.")
                dir_to_save = self.working_directory
                filepath = path.join(dir_to_save, filename)
                with open(file=filepath, mode="wb") as fobj:
                    pickle.dump(self, file=fobj)
            else:
                raise exc
        return filepath
    
    @staticmethod
    def load_snapshot_file(filepath: str) -> ExecutionEntity:
        """
        Loads a ExecutionEntity instance from a snapshot file.

        This method deserializes the state of a ExecutionEntity instance from the given
        snapshot file. If the deserialized object is not an instance of ExecutionEntity,
        it raises a TypeError.

        Args:
            filepath (str): The path to the snapshot file from which to load the state.

        Returns:
            ExecutionEntity: The deserialized ExecutionEntity instance.

        Raises:
            TypeError: If the deserialized object is not an instance of ExecutionEntity.
        """
        with open(file=filepath, mode="rb") as fobj:
            unit = pickle.load(file=fobj)
            
            if not isinstance(unit, __class__):
                raise TypeError(f"The loaded object is not an instance of {__class__.__name__}.")
            
            database_session = unit.create_database_session()
            database_session.close()
            return unit

    def reload(self):
        """Reload task. Placeholder ONLY!"""
        if self.status in StatusCode.unexecuted_statuses:
            _LOGGER.info(f"{self.identifier} hasn't been executed yet, cannot reload!")
            return
        _LOGGER.info(f"Reloading {self.identifier}...")

class WorkFlow(ExecutionEntity):
    """Rrepresents a procedural WorkFlow primarily consists of a number of WorkUnit instances, 
    including data mappers, self-inspection flags, to ensure its operation.

    Each WorkFlow Instance will perform self-inspection during initialization. If the self-inspection fails, it cannot be executed.
    
    If an error occurs during the execution of a single WorkUnit instance, the workflow will record the error, 
    jump out of the WorkUnit instance, and execute the next one until the execution of the entire workflow is completed.
    
    Attributes:
        api_key (str): The value of api_key of workflow instance, which is set to be WORKFLOW_API_KEY.
        error_msg_list (list): A list of error logs collected from the initialization and execution of the workflow. This is a class attribute.
        control_workunit (ControlWorkUnit): The outer ControlWorkUnit instance hosting the current workflow.
        workunits (List[WorkFlow]): A list of WorkUnit instance carried by this workflow.
        intermediate_data_mapper (dict): Runtime intermediate data for the current layer of workflow.
        inherited_data_mapper (dict): Data inherited from parent workflow.
        nested_control_unit_return_keys (list): The return keys of nested ControlWorkUnit instances (e.g. loop)
                                                Keys stored in it will not be inherited to prevent the `inherited_data_mapper` from being toooo large.
        locator (List[int]): A list of the hierarchical indexes to locate the current workflow.
                                  If the current workflow is at the outermost level, this value is an empty list.
        execution_status (int): A flag indicating the execution status of the WorkFlow instance. Using value from `class ExecutionStatus`.
        Other inherited attributes...
    """
    api_key: str
    error_msg_list: list     # Error information. This is a class attribute.
    
    control_workunit: ControlWorkUnit   # The outer ControlWorkUnit instance hosting the current workflow.
    workunits: List[WorkUnit]           # A list of WorkUnit instance carried by this workflow.
    
    # Intermediate Data.
    intermediate_data_mapper: dict  # Runtime intermediate data for the current layer of workflow.
    inherited_data_mapper: dict     # Data inherited from parent workflow.
    nested_control_unit_return_keys: list   # The return keys of nested ControlWorkUnit instances (e.g. loop)

    def __init__(self, working_directory: str = BASE_DIRECTORY, 
                debug: bool = False, 
                locator: List[int] = list(), 
                intermediate_data_mapper: Dict[str, Any] = dict(), 
                inherited_data_mapper: Dict[str, Any] = dict(),
                control_workunit: ControlWorkUnit = None,
                sqlite_filepath: str = str(),
                ) -> None:
        """Initialize an instance.
        
        Args:
            debug (bool, optional): Indicates whether to run in debug mode.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (List[int]): A list of the hierarchical indexes to locate the current workflow.
                                If the current workflow is at the outermost level, this value is an empty list.
            intermediate_data_mapper (Dict[str, Any], optional): The initial value of intermediate_data_mapper.
            inherited_data_mapper (Dict[str, Any], optional): Data Mapper containing data to be inherited.
            control_workunit (ControlWorkUnit): The outer ControlWorkUnit instance hosting the current workflow.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
        """
        super().__init__(working_directory=working_directory, debug=debug, locator=locator, sqlite_filepath=sqlite_filepath)

        # Assigned values.
        self.intermediate_data_mapper = intermediate_data_mapper
        self.inherited_data_mapper = inherited_data_mapper
        _LOGGER.debug(f"inherited_data_mapper: {self.inherited_data_mapper}")
        self.control_workunit = control_workunit

        # Default values.
        self.api_key = WORKFLOW_API_KEY
        self.error_msg_list = list()
        self.workunits = list()
        self.nested_control_unit_return_keys = list()
        self.child_execution_entities = self.workunits
        return
    
    @classmethod
    def from_list(cls, unit_dict_list: List[dict], working_directory: str = BASE_DIRECTORY, 
                debug: bool = False, locator: List[int] = list(),
                intermediate_data_mapper: Dict[str, Any] = dict(), 
                inherited_data_mapper: Dict[str, Any] = dict(), 
                control_workunit: ControlWorkUnit = None, 
                sqlite_filepath: str = str(),
                ) -> WorkFlow:
        """Initialize an instance of WorkFlow from a given dictionary instance.
        
        Args:
            unit_dict_list (List[dict]): A list containing the WorkUnit configuration dicts.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            debug (bool, optional): Indicates whether to run in debug mode.
            locator (List[int]): A list of the hierarchical indexes to locate the current workflow.
                                If the current workflow is at the outermost level, locator is an empty list.
            intermediate_data_mapper (Dict[str, Any], optional): The initial value of intermediate_data_mapper.
            inherited_data_mapper (Dict[str, Any], optional): Data Mapper containing data to be inherited.
            control_workunit (ControlWorkUnit): The outer ControlWorkUnit instance hosting the current workflow.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
        
        Returns:
            An instance of WorkFlow.
        """
        flow = cls(working_directory, debug, locator, intermediate_data_mapper, inherited_data_mapper, control_workunit, sqlite_filepath)
        flow.status = StatusCode.INITIALIZING
        if (unit_dict_list == None):
            no_workunit_err_msg = "Initializing WorkFlow with no workunits. Workunits are expected."
            flow.error_msg_list.append(no_workunit_err_msg)
            flow.status = StatusCode.FAILED_INITIALIZATION
        else:
            for unit_index, unit_dict in enumerate(unit_dict_list):
                locator_to_pass = locator.copy()
                locator_to_pass.append(unit_index)
                flow.add_unit(unit_dict, locator_to_pass)
                continue
        
        # If the flow status keeps INITIALIZING, then mark it as READY_TO_START.
        if flow.status == StatusCode.INITIALIZING:
            flow.status = StatusCode.READY_TO_START
        return flow
    
    @classmethod
    def from_json_string(cls, json_str: str, working_directory: str = BASE_DIRECTORY, debug: bool = False) -> WorkFlow:
        """Initialize an instance of WorkFlow from a serialized json string.
        
        Args:
            json_str (str): The serialized json string.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        unit_dict_list = loads(json_str)
        return WorkFlow.from_list(unit_dict_list, working_directory=working_directory, debug=debug)

    @classmethod
    def from_json_file_object(cls, json_fobj: TextIOWrapper, working_directory: str = BASE_DIRECTORY, debug: bool = False) -> WorkFlow:
        """Initialize an instance of WorkFlow from a json file object.
        
        Args:
            json_fobj (TextIOWrapper): The instance of a json file containing workflow configuration information.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        unit_dict_list = load(json_fobj)
        return WorkFlow.from_list(unit_dict_list, working_directory=working_directory, debug=debug)

    @classmethod
    def from_json_filepath(cls, json_filepath: str, working_directory: str = BASE_DIRECTORY, debug: bool = False) -> WorkFlow:
        """Initialize an instance of WorkFlow from a json filepath.
        
        Args:
            json_filepath (str): The path to a json file containing workflow configuration information.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        with open(json_filepath) as fobj:
            return WorkFlow.from_json_file_object(fobj, working_directory=working_directory, debug=debug)
    
    def add_unit(self, unit_dict: dict, locator_to_pass: List[int] = list()) -> None:
        """Adds a new unit to the workflow.

        This method creates a new instance of WorkUnit from a provided dictionary and adds it to 
        the workflow. It initializes the WorkUnit with the given configuration `unit_dict`, 
        links it to the current workflow instance, and sets the debug mode according to the 
        provided argument. Any error messages encountered during the initialization of the WorkUnit 
        are appended to the workflow's error message list.

        Args:
            unit_dict (dict): A dictionary containing the configuration for the unit to be added.
                            This dictionary should include keys for API mapping and execution 
                            parameters.
            locator_to_pass (List[int]): A list of the hierarchical indexes to locate the new workunit.
        
        Returns:
            The WorkUnit instance added.
        """
        api_key = unit_dict.get(WORKUNIT_API_NAME_KEY)
        if (api_key in CONTROL_API_KEYS):
            _LOGGER.info(f"Parsing Control WorkUnit {api_key}.")
            if api_key == LOOP_API_KEY:
                # Setup a placeholder in `workflow.intermediate_data_mapper`.
                workunit = LoopWorkUnit.from_dict(unit_dict=unit_dict, workflow=self, working_directory=self.working_directory, locator=locator_to_pass, sqlite_filepath=self.sqlite_filepath, debug=self.debug)
            elif api_key == CLUSTER_BATCH_API_KEY:
                workunit = ClusterBatchWorkUnit.from_dict(unit_dict=unit_dict, workflow=self, working_directory=self.working_directory, locator=locator_to_pass, sqlite_filepath=self.sqlite_filepath, debug=self.debug)
            elif api_key == GENERAL_API_KEY:
                error_msg = "GeneralWorkUnit should not be configured inside a workflow."
                _LOGGER.error(error_msg)
                self.error_msg_list.append(error_msg)
        else:
            _LOGGER.info(f"Parsing Science WorkUnit {api_key}.")
            workunit = WorkUnit.from_dict(unit_dict=unit_dict, workflow=self, working_directory=self.working_directory, locator=locator_to_pass, sqlite_filepath=self.sqlite_filepath, debug=self.debug)

        # Setup a placeholder in `workflow.intermediate_data_mapper`.
        self.intermediate_data_mapper[VarFormatter.unformat_variable(workunit.return_key)] = Placeholder.PLACEHOLDER_STR_VALUE

        # Add new workunit to the workflow.
        self.workunits.append(workunit)
        if workunit.status == StatusCode.FAILED_INITIALIZATION:
            self.status = StatusCode.FAILED_INITIALIZATION
            self.error_msg_list += workunit.error_msg_list
        return

    def reload(self, unit_dict_list: List[dict]):
        """
        Reloads the WorkFlow instance with updated configuration.

        This method is used when a user wants to reload the entire task and checks 
        if there are any changes in the parameter values. To perform this, we reload each
        WorkUnit instance in this WorkFlow instance and update its intermediate_data_mapper.

        Args:
            unit_dict_list (List[dict]): A list containing the WorkUnit configuration dicts.

        Notes:
            - If a WorkUnit instance is marked as `SUSPECIOUS_UPDATE` after reloading, 
                its `store_as` value will be assigned as a placeholder value in the intermediate_data_mapper.
        """
        super().reload()
        has_suspecious_updates = False
        has_unresolved_error_or_pause = False
        for i, unit_dict in enumerate(unit_dict_list):
            workunit = self.workunits[i]
            store_as_key_backup = workunit.return_key
            workunit.reload(unit_dict)
            if workunit.status in StatusCode.unexecutable_statuses:
                self.status = StatusCode.FAILED_INITIALIZATION
            elif workunit.status == StatusCode.SUSPECIOUS_UPDATES:
                has_suspecious_updates = True
                store_as_key = workunit.return_key
                store_as_value = workunit.return_value
                del self.intermediate_data_mapper[store_as_key_backup]
                self.intermediate_data_mapper[store_as_key] = Placeholder.assign_placeholder_value(store_as_value)
            elif workunit.status in StatusCode.error_or_pause_statuses:
                has_unresolved_error_or_pause = True
            continue
        if self.status not in StatusCode.unexecutable_statuses:
            if has_unresolved_error_or_pause:
                return
            elif has_suspecious_updates:
                self.status = StatusCode.SUSPECIOUS_UPDATES
            else:
                self.status = StatusCode.EXIT_OK

    def execute(self):
        """Executes all units in the workflow sequentially.

        Before execution, this method checks if the workflow has passed its self-inspection. 
        If the self-inspection has not been passed, it logs an error and raises a ValueError, 
        indicating that the workflow is not ready for execution. If the self-inspection is 
        passed, it proceeds to execute each work unit in the order they were added to the 
        workflow. 

        During the execution of each unit, the method captures and stores the return values 
        and any error messages generated. The return values are stored in the 
        `workunits_data_mapper` dictionary, mapping the return keys to their corresponding 
        return values. Any error messages generated during the execution of work units are 
        appended to the workflow's error message list.

        This process design of "logging errors" and holding off on "raising errors"
        is conducive to keeping the workflow from being interrupted by non-fatal errors,
        so that the workflow can continue to run until completion as long as there are no fatal errors.

        Note: Chances exist that the `intermediate_data_mapper` may have `None: None` key-value pairs,
        but that's fine: A `None: None` key-value pair wouldn't make any difference. Thus, we don't need
        to add an if-else to filter it out.

        Raises:
            ValueError: If the workflow has not passed the self-inspection.
        """
        if not self.is_executable:
            _LOGGER.error(f"{self.identifier}: {FAILED_INITIALIZATION_ERROR_MSG}")
            return dict()
        else:
            if not self.control_workunit:    # If the current workflow is at the outermost level, log the message; otherwise skip the duplicated message.
                _LOGGER.info(SUCCESS_INITIALIZATION_MSG)
            self.status = StatusCode.RUNNING
            for workunit in self.workunits:
                workunit: WorkUnit
                _LOGGER.info(f"Executing WorkUnit {workunit.identifier} ...")
                return_key, return_value = workunit.execute()
                self.intermediate_data_mapper[VarFormatter.unformat_variable(return_key)] = return_value
                if (workunit.status in StatusCode.error_excluding_pause_statuses):
                    self.error_msg_list += workunit.error_msg_list
                    self.status = StatusCode.EXIT_WITH_ERROR_IN_INNER_UNITS
                    continue
                elif (workunit.status in StatusCode.pause_excluding_error_statuses):
                    if (self.status not in StatusCode.error_including_pause_statuses):
                        self.status = StatusCode.RUNNING_WITH_PAUSE_IN_INNER_UNITS
                    break
                elif (workunit.status in StatusCode.pause_including_error_statuses):
                    self.status = StatusCode.EXIT_WITH_ERROR_AND_PAUSE
                    break
            if self.status == StatusCode.RUNNING:
                self.status = StatusCode.EXIT_OK
            return self.intermediate_data_mapper
        
    def execute_cluster(self):
        """Submit and execute the WorkFlow instance in a computational cluster."""
        pass

class WorkUnit(ExecutionEntity):
    """Represents a single unit in a procedural workflow.

    - This class encapsulates the information and functionality for initializing, interpreting and executing
    a single unit of work in a procedural workflow. 
    - It handles API mappings, argument processing, and execution
    of a specified function or method.

    Attributes:
        api_key (str): Key used to identify and map the appropriate science API function.
                         Read from the `api` field in the configuration.
        args_dict_input (dict): The raw, unprocessed arguments dictionary as specified
                                  in each unit of the JSON configuration.
        workflow (WorkFlow): Reference to the parent workflow object that this unit is part of.
        api (object): The actual science API function object mapped using `__api_key`.
        return_key (str): The key under which the return value of the API call will be stored.
                        Read from the `store_as` field in the configuration.
        return_value: The value returned by the execution of the API function.
        args_dict_to_pass (dict): Processed arguments dictionary to be passed to the API function.
        args_dict_to_pass_backup (dict): For reload use only. A backup of `args_dict_to_pass` to compare before execution.
        error_info_list (list): List of error messages encountered during processing.
        params_to_assign_at_execution (list): List of parameters that need to be assigned values at execution time.
        architecture (str): The workunits and workflows to be performed in the task, and their order and hierarchy of operation.
        Other inherited attributes...
    """
    api_key: str
    args_dict_input: dict
    architecture: str

    workflow: WorkFlow
    api: Callable
    return_key: str
    return_value: Any
    args_dict_to_pass: dict
    args_dict_to_pass_backup: dict
    error_msg_list: list
    params_to_assign_at_execution: list

    def __init__(self, unit_dict: dict = dict(), workflow: WorkFlow = None, working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), sqlite_filepath: str = str(), debug: bool = False):
        """Initializes a WorkUnit instance.

        Args:
            unit_dict (dict): A dictionary derived from a JSON configuration, containing keys
                              for API mapping and execution.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (list, optional): A list of the hierarchical indexes to locate the current workunit.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        super().__init__(working_directory=working_directory, debug=debug, locator=locator, sqlite_filepath=sqlite_filepath)
        # Initialize attributes in the constructor
        # to prevent them from being considered as class attributes.
        self.architecture = __class__.read_architecture(unit_dict=unit_dict)
        self.api_key = unit_dict.get(WORKUNIT_API_NAME_KEY, str())
        self.return_key = unit_dict.get(WORKUNIT_RETURN_VALUE_KEY, None)
        self.args_dict_input = unit_dict.get(WORKUNIT_ARGUMENT_LIST_KEY, dict())
        self.args_dict_to_pass = dict()
        self.error_msg_list = list()
        self.params_to_assign_at_execution = list()
        self.return_value = None

        if workflow == None:
            self.workflow = WorkFlow()
        else:
            self.workflow = workflow
        return
    
    @staticmethod
    def read_architecture(unit_dict: dict, indent: int = 0):
        """
        Recursively reads a nested unit_dict dictionary and extracts 'api' values.

        Args:
            schema (dict): The nested dictionary schema.
            indent (int): The current indentation level (number of spaces).

        Returns:
            str: A formatted string representing the 'api' values at each level.
        """
        result = ""

        # Check if 'api' key exists in the current level
        if WORKUNIT_API_NAME_KEY in unit_dict:
            result += ' ' * indent + f"- {unit_dict[WORKUNIT_API_NAME_KEY]}\n"
        
        # If 'args' key exists and it has a 'workunits' key, process each item in 'workunits'
        if (WORKUNIT_ARGUMENT_LIST_KEY in unit_dict) and (CONTROL_BODY_WORKUNITS_LABEL in unit_dict[WORKUNIT_ARGUMENT_LIST_KEY]):
            for workunit in unit_dict[WORKUNIT_ARGUMENT_LIST_KEY][CONTROL_BODY_WORKUNITS_LABEL]:
                result += __class__.read_architecture(workunit, indent + 4)

        return result

    @classmethod
    def from_dict(cls, unit_dict: dict, workflow: WorkFlow = None, working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), sqlite_filepath: str = str(), debug: bool = False) -> WorkUnit:
        """Initializes an instance of WorkUnit from a given dictionary.

        The dictionary should contain keys such as `api`, `store_as`, and `args`
        specifying the API to be called, the return key for storing the result, and the arguments
        for the API call, respectively.

        Args:
            unit_dict (dict): A dictionary derived from a JSON configuration, containing keys
                              for API mapping and execution.
            workflow (WorkFlow, optional): The workflow instance to which this unit belongs.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (list, optional): A list of the hierarchical indexes to locate the current workunit.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            WorkUnit: An initialized instance of WorkUnit.
        
        Raises:
            KeyError: If the API key cannot be mapped to any API.
        """
        # Initialize the class.
        unit = cls(unit_dict=unit_dict, workflow=workflow, working_directory=working_directory, locator=locator, sqlite_filepath=sqlite_filepath, debug=debug)
        unit.status = StatusCode.INITIALIZING

        # Map the API.
        unit.api = SCIENCE_API_MAPPER.get(unit.api_key, None)
        if unit.api is None:    # Inspect API Mapping.
            error_msg = f"Input API Key `{unit.api_key}` cannot be mapped to any API."
            _LOGGER.error(error_msg)
            unit.status = StatusCode.FAILED_INITIALIZATION
            if unit.debug:
                raise KeyError(unit.api_key)
            else:
                unit.error_msg_list.append(error_msg)
                return unit
        
        _LOGGER.debug(f"The API of '{unit.identifier}' is '{unit.api}'.")
        # Self Inspection.
        unit.self_inspection_and_args_reassembly()
        unit.check_initialization_status()
        return unit

    def reload(self, unit_dict: dict, api: Callable = None) -> WorkUnit:
        """
        Reloads the WorkUnit with updated configuration.

        This method is used when a user wants to reload the entire task and checks 
        if there are any changes in the parameter values. If the WorkUnit's status 
        indicates an error or pause, it updates the WorkUnit's properties with the 
        new configuration. Otherwise, it compares the current and new configurations 
        to check for updates.

        Args:
            unit_dict (dict): The new configuration dictionary for the WorkUnit.

        Notes:
            - Compares `args_dict_input`, `args_dict_to_pass`, and `params_to_assign_at_execution`
              between the current WorkUnit and a virtual unit created from the new configuration.
            - Also compares `return_key` values.
            - If any of these values are different, sets the status to `StatusCode.SUSPECIOUS_UPDATES`.

        Returns:
            The created virtual WorkUnit instance for further comparison.
        """
        super().reload()
        # Create a virtual unit with the new configuration for comparison.
        self.args_dict_to_pass_backup = self.args_dict_to_pass.copy()
        
        virtual_unit = WorkUnit(unit_dict=unit_dict)
        virtual_unit.status = StatusCode.INITIALIZING
        virtual_unit.workflow = self.workflow
        if (api is None):
            virtual_unit.api = SCIENCE_API_MAPPER.get(virtual_unit.api_key)
        else:
            virtual_unit.api = api
        virtual_unit.self_inspection_and_args_reassembly()
        virtual_unit.check_initialization_status()
        updated_params = virtual_unit.update_args_value_from_data_mapper()

        if (virtual_unit.status in StatusCode.unexecutable_statuses):
            self.status = StatusCode.FAILED_INITIALIZATION
            self.check_initialization_status()

        # Compare the attributes for updates.
        # Check if any of the attributes (args_dict_input, args_dict_to_pass, params_to_assign_at_execution,
        # or return_key) are different in the virtual unit.
        # If there are differences, update the attributes and then update status to indicate that updates are detected.
        elif (updated_params or
            self.return_key != virtual_unit.return_key or 
            self.args_dict_input != virtual_unit.args_dict_input or 
            self.args_dict_to_pass != virtual_unit.args_dict_to_pass or 
            self.params_to_assign_at_execution != virtual_unit.params_to_assign_at_execution):
            
            self.status = StatusCode.SUSPECIOUS_UPDATES
            self.args_dict_input = virtual_unit.args_dict_input
            self.args_dict_to_pass = virtual_unit.args_dict_to_pass
            self.params_to_assign_at_execution = virtual_unit.params_to_assign_at_execution
            self.return_key = virtual_unit.return_key
        
        return virtual_unit

    def execute(self) -> Tuple[str, Any]:
        """Executes the workunit.

        This method attempts to execute the API function associated with this work unit using the
        arguments prepared during initialization. The result of the execution is stored and mapped
        with its return key. If an exception occurs during execution, the behavior depends on the
        debug mode: if debug is True, the exception is raised immediately; otherwise, it logs each
        error and record errors in the error list.

        Returns:
            tuple: A tuple containing the `return_key` and the `return_value` from the API execution.
                If an error occurs and debug is False, `return_value` will be None.

        Raises:
            Exception: Any exception raised by the API function if debug mode is enabled.
        """
        try:
            self.update_args_value_from_data_mapper()
            _LOGGER.debug(f"{self.identifier} with Status {self.status}.")
            # if (self.status not in StatusCode.unexecuted_statuses and self.args_dict_to_pass == self.args_dict_to_pass_backup):
            #     if (self.status == StatusCode.SUSPECIOUS_UPDATES):
            #         self.status = StatusCode.EXIT_OK
            #     elif (self.status in StatusCode.error_or_pause_statuses):
            #         # If the old status is error or pause 
            #         # and nothing changed to `args_dict_to_pass`, no need to run it again.
            #         return self.return_key, self.return_value

            if (self.status in StatusCode.skippable_statuses):
                _LOGGER.info(f"WorkUnit {self.identifier} is skipped!")
                return self.return_key, self.return_value
                
            self.status = StatusCode.RUNNING

            if self.return_key:
                self.return_value = self.api(**self.args_dict_to_pass)
            else:   # If no `return_key`, then `return_value` is None.
                self.api(**self.args_dict_to_pass)
                self.return_value = None
            self.status = StatusCode.EXIT_OK
            _LOGGER.debug(f"{self.identifier}(return_key={self.return_key}, return_value={self.return_value})")
            return self.return_key, self.return_value
        except Exception as exc:
            _LOGGER.error(f'Catching Error when executing `{self.identifier}`:')
            self.status = StatusCode.EXIT_WITH_ERROR
            if (self.debug):
                raise exc
            else:
                _LOGGER.error(exc)
                self.error_msg_list.append(str(exc))
                return self.return_key, None

    def update_args_value_from_data_mapper(self, commit_update: bool = True) -> list:
        """Update the `args_dict_to_pass` by mapping data from `workflow.intermediate_data_mapper` or `workflow.inherited_data_mapper`.
        If a value is mapped in both `inherited_data_mapper` and `intermediate_data_mapper`,
        then the latter one will overwrite the former one.

        Args:
            commit_update (bool, optional): Whether to update discovered parameter value changes to `args_dict_to_pass`?

        Returns:
            A list of updated params.
        """
        updated_params = list()
        for param in self.params_to_assign_at_execution:
            for mapper in [self.workflow.intermediate_data_mapper, self.workflow.inherited_data_mapper]:
                if (mapped_datum := mapper.get(VarFormatter.unformat_variable(self.args_dict_input[param]))):
                    existing_arg_in_mapper = self.args_dict_to_pass.get(param)
                    if (Placeholder.is_none_or_placeholder(existing_arg_in_mapper) or (existing_arg_in_mapper != mapped_datum)):
                        if (commit_update):
                            if (Placeholder.is_none_or_placeholder(mapped_datum)):
                                continue
                            else:
                                self.args_dict_to_pass[param] = mapped_datum
                        updated_params.append(param)
                        break
            continue
        return updated_params
    
    @staticmethod
    def __inspect_argument_type(arg_value: Any, annotation: Union[Type, Any]) -> bool:
        """
        Inspects whether the given argument value matches the specified type annotation.

        This method supports both basic types and typing types (e.g., `List`, `Dict` from `typing` module).
        For typing types, it checks if the argument value matches the origin type of the annotation
        (e.g., list for List[int], dict for Dict[str, int]).
        Note: Generic elements in typing types, like `int` in `List[int]`, `str`
        and `int` in `Dict[str, int]`, are not inspected, as JSON files do not support storing
        such types directly (Zhong: TODO later).

        Due to some unknown reason, some annotations may not be parse correctly. In such cases, we will set
        the return value as True with a WARNING logged so as to ensure the smooth operation of workflow.

        Args:
            arg_value (Any): The value of the argument to be inspected.
            annotation (Union[Type, Any]): The type annotation against which the argument value
                                           is to be checked. This can be a basic type like `int`,
                                           `str`, or a typing type like `List[str]`.

        Returns:
            bool: True if the argument value matches the type annotation, False otherwise.
        """
        if (annotation == Any):
            return True
        elif hasattr(annotation, '__origin__'):   # Checks if param.annotation is a typing type.
            if (annotation.__origin__ is Union):   # Handle Union Type.
                valid_types = get_args(annotation)
                _LOGGER.debug(f"Type from Union are {valid_types}.")
                subscripted_generics = list()
                for typ in valid_types:
                    try:
                        if isinstance(arg_value, typ):
                            return True 
                    except (TypeError):
                        subscripted_generics.append(typ)
                        continue
                if (subscripted_generics):
                    _LOGGER.warning(f"Please make sure that `{arg_value}` is an instance of {' or '.join(subscripted_generics)}.")
                    return True
                else:
                    return False
            else:
                # Check the origin of typing types (e.g., list for List[int])
                return isinstance(arg_value, annotation.__origin__)
        else:
            if (not isinstance(annotation, type)):
                # TODO Sometimes, the annotation may be parsed as a string value,
                # so we use pydoc.locate to cast it from string to type if it happens.
                if (isinstance(annotation, str)):
                    index_of_bracket = annotation.find("[")
                    matched_type = annotation[:index_of_bracket].lower() if index_of_bracket != -1 else annotation.lower()
                    from pydoc import locate
                    matched_type = locate(matched_type)

                    if matched_type:
                        annotation = matched_type
                    else:
                        annotation = eval(annotation)
                else:
                    _LOGGER.warning(f"Unable to determine if `{arg_value}` is an instance of {annotation}. Please ensure it by yourself")
                    return True

            # Inspect basic types.
            return isinstance(arg_value, annotation)

    def __load_data_mapper_value(self, param_name: str, arg_value: str) -> Tuple[bool, Any]:
        """Load data from `intermediate_data_mapper` or `inherited_data_mapper`.
        
        Returns:
            A tuple containing a boolean value of "whether to skip type inspection" and loaded `arg_value`.
        """
        unformatted_arg_value = VarFormatter.unformat_variable(arg_value)
        is_key_mapped = False   # Indicate if the value is mapped.
        is_value_mapped_none = True # Indicate if the value in the mapper is None at present.
        for mapper in [self.workflow.intermediate_data_mapper, self.workflow.inherited_data_mapper]:    # Listed in order of priority.
            if (mapper and unformatted_arg_value in mapper):    # If mapped, leave for later.
                is_key_mapped = True
                self.params_to_assign_at_execution.append(param_name)
                self.args_dict_to_pass[param_name] = unformatted_arg_value
                mapped_value = mapper.get(unformatted_arg_value)
                _LOGGER.debug(f"Checkpoint 1: The mapped value is `{mapped_value}`.")
                if (mapped_value != None):   # Check if there is a non-none-or-empty value in Mapper.
                    _LOGGER.debug(f"Checkpoint 2: The mapped value `{mapped_value}` is not NoneValue.")
                    if (isinstance(mapped_value, type(Placeholder.PLACEHOLDER_STR_VALUE)) # Placeholder values need to be excluded.
                        and mapped_value == Placeholder.PLACEHOLDER_STR_VALUE):
                        _LOGGER.debug(f"Checkpoint 3: The mapped value `{mapped_value}` is a placeholder value.")
                        break
                    is_value_mapped_none = False
                    arg_value = mapped_value    # Assign mapped value to the argument.
                if is_key_mapped:   # Once key is mapped, no longer mapping other mappers.
                    break
        if (is_key_mapped and is_value_mapped_none):    # If the key is mapped but the value is still None, skip inspection.
            return True, arg_value
        elif (not is_key_mapped):
            error_msg = f"The param `{param_name}` in API `{self.api_key}` hasn't received expected variable {VarFormatter.format_as_variable(arg_value)}. Check your configuration and there's maybe some typos!"
            _LOGGER.error(error_msg)
            self.error_msg_list.append(error_msg)
            return True, arg_value
        _LOGGER.debug(f"{self.api_key}, param_name: {param_name}, arg_value: {unformatted_arg_value}, Value Changed.")
        return False, arg_value

    def self_inspection_and_args_reassembly(self) -> None:
        """Performs self-inspection and reassembles the arguments for the API call.

        This method is called every time a new instance is created. It performs a self-inspection
        to validate the `args` dictionary specified in the JSON configuration. It checks for the
        existence of required arguments and their types. If any errors are found, they are recorded
        in `error_info_list`, and the `status` flag is set to `FAILED_INITIALIZATION`.

        The arguments that are to be mapped from the mapper are reassembled into `args_dict_to_pass`
        for the API call. Such arguments, if they already have been assigned with value in the mapper,
        they would have to undergo the inspection; otherwise (if is None or Empty value), they could skip the inspection.
        """
        unit_error_list = list()
        
        # Inspect Arguments.
        api_parameters = signature(self.api).parameters
        for param_name in api_parameters:
            param: Parameter = api_parameters[param_name]

            # If `*args` or `**kwargs` parameters exist (which is always at the very last),
            # add all the arguments to the `args_dict_to_pass`. Then, break the loop.
            # APIs with both `*args` and `**kwargs` are not and will not be supported.
            if (param_name in ['arg', 'args', 'kwarg', 'kwargs']):
                for key, value in self.args_dict_input.items():
                    if key not in self.args_dict_to_pass.keys():
                        if VarFormatter.is_variable_formatted(value):   # Check if the arg_value is a formatted as a value name.
                            _, arg_value = self.__load_data_mapper_value(key, value)  # kwargs params always skip type inspection.
                            self.args_dict_to_pass[key] = arg_value
                        else:
                            self.args_dict_to_pass[key] = value
                            continue
                break

            # A param with empty default value is a required parameter.
            is_required = False
            if param.default == _empty:
                is_required = True

            # Inspect argument existence.
            arg_value = self.args_dict_input.get(param.name, None)
            if not arg_value:   # Missing the argument.
                if is_required: # Record error if it's required.
                    error_msg = f'Missing required argument `{param.name}` in API `{self.api_key}`.'
                    unit_error_list.append(error_msg)
                    _LOGGER.error(error_msg)
                    self.status = StatusCode.FAILED_INITIALIZATION
                else:   # Skip if it's optional.
                    pass
                continue

            # Load data from `intermediate_data_mapper` or `inherited_data_mapper`.
            _LOGGER.debug(f"{self.api_key}, param_name: {param.name}, arg_value: {arg_value}")
            skip_type_inspection = False
            if VarFormatter.is_variable_formatted(arg_value):   # Check if the arg_value is a formatted as a value name.
                skip_type_inspection, arg_value = self.__load_data_mapper_value(param_name, arg_value)
            
            if skip_type_inspection:    # Skip type inspection if `skip_type_inspection` is True.
                continue

            # Inspect argument type.
            if (param.annotation != _empty):
                annotation = param.annotation
                pass_type_inspection = WorkUnit.__inspect_argument_type(arg_value, annotation)
                if not pass_type_inspection:
                    error_msg = f'Receiving argument `{param.name}` with value `{arg_value}` in unexpected type {type(arg_value)} while {annotation} is expected (when initializing `{self.api_key}`).'
                    _LOGGER.error(error_msg)
                    unit_error_list.append(error_msg)
                    self.status = StatusCode.FAILED_INITIALIZATION

            # Inspection passed, now loading arguments.
            self.args_dict_to_pass[param.name] = arg_value
            continue
        
        self.error_msg_list += unit_error_list
        return
    
    def check_initialization_status(self) -> None:
        """Checks the result of self-inspection of the instance.

        This method checks if the initialization fails (`self.status`). 
        If the initialization fails and debugging is enabled, it logs each error message and raises a 
        ValueError. If debugging is not enabled, nothing happens.
        If the initialization succeeds, then the status of the instance is set to READY_TO_START.

        Raises:
            ValueError: If initialization fails and debugging is enabled, a ValueError is raised
                        with the input argument dictionary.
        """
        if self.status == StatusCode.FAILED_INITIALIZATION:
            if self.debug:
                raise ValueError(self.args_dict_input)
            else:
                return
        else:
            self.status = StatusCode.READY_TO_START
            return

    def generate_data_mapper_to_inherite(self) -> Dict[str, Any]:
        """Generate data mappers for sub-workflow inheritance.
        
        Returns:
            Data mapper to inherite (dict)
        """
        data_mapper_to_inherite = dict()
        if (inherited_data_mapper:=self.workflow.inherited_data_mapper):
            # Data Mapper inherited from parent workflow(s) should be inherited.
            data_mapper_to_inherite.update(inherited_data_mapper.copy())
        for key, value in self.workflow.intermediate_data_mapper.items():
            if (self.workflow and key in self.workflow.nested_control_unit_return_keys):
                # Keys stored in `nested_control_unit_return_keys` should not be inherited
                continue
            data_mapper_to_inherite[key] = value
        return data_mapper_to_inherite

#endregion

#region Non-entry Control entities.

class ControlWorkUnit(WorkUnit):
    """This class is an abstract class defining shared attributes and methods of Control WorkUnits (e.g. loop, general, etc.)
    This class is derived from WorkUnit class.
    """

    def __init__(self, unit_dict: dict = dict(), workflow: WorkFlow = None, 
                working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), 
                sqlite_filepath: str = str(), debug: bool = False):
        super().__init__(unit_dict=unit_dict, workflow=workflow, working_directory=working_directory, 
                    locator=locator, sqlite_filepath=sqlite_filepath, debug=debug)
        if (self.workflow):
            self.workflow.nested_control_unit_return_keys.append(self.return_key)

    def reload(self, unit_dict: dict, api: Callable) -> ControlWorkUnit:
        """
        Reloads the Control WorkUnit instance with updated configuration.

        This method is used when a user wants to reload the entire task and checks 
        if there are any changes in the parameter values. If the WorkUnit's status 
        indicates an error or pause, it updates the WorkUnit's properties with the 
        new configuration. Otherwise, it compares the current and new configurations 
        to check for updates.

        Args:
            unit_dict (dict): The new configuration dictionary for the WorkUnit.
            api (Callable): The Callable object to perform self-inspection.

        Notes:
            - Compares `args_dict_input`, `args_dict_to_pass`, and `params_to_assign_at_execution`
              between the current WorkUnit and a virtual unit created from the new configuration.
            - Also compares `return_key` values.
            - If any of these values are different, sets the status to `StatusCode.SUSPECIOUS_UPDATES`.

        Returns:
            The created virtual WorkUnit instance for further comparison.
        """
        virtual_unit = super().reload(unit_dict=unit_dict, api=api)
        return virtual_unit

    def encode_layer_index(self, execution_index: int = None) -> str:
        """The index of current layer to be added to `locator`.
        
        Args:
            execution_index (int, optional): Indicates the execution index of the target sub-workflow instance in the current control unit, 
                                            such as the index of loops in the `LoopWorkUnit`,
                                            and/or the index of parallel tasks in the `ClusterBatchWorkUnit`.
        
        Returns:
            The encoded index of current layer. Default `self.api_key`
        """
        splitter = "_"
        if execution_index is None:
            return self.api_key
        else:
            return f"{self.api_key}{splitter}{execution_index}"
    
    def decode_layer_index(self, layer_index: str) -> int:
        """Decode the index of current layer to obtain the execution index.
        
        Args:
            layer_index (str): The index of current layer in `locator`.

        Return:
            The execution index obtained. Default 0.
        """
        api_key_length = len(self.api_key)
        splitter = "_"
        if layer_index[:api_key_length] == self.api_key and layer_index[api_key_length] == splitter:
            if (splitter in layer_index):
                index_str = layer_index.split(splitter)[-1]
                index = int(index_str)
                return index
        else:
            return 0
    
    def locate(self, relative_locator: list) -> ExecutionEntity:
        """Recursive locating of ExecutionEntity instance based on their locators.

        Note: This method will return None for an unexecuted ExecutionEntity instance that is inside a ControlWorkUnit instance other than GeneralWorkUnit.
        
        Args:
            relative_locator (list): A locator that locates the target instance starting from the current instance.
                                    For example, if the locator of target instance is [4, loop_0, 1, 2] and the locator of current instance is [4, loop_0],
                                    then the value of `relative_locator` should be [1, 2].

        Returns:
            The located ExecutionEntity instance. None if nothing is located.
        """
        if len(relative_locator) == 0:
            return self
        else:
            try:
                located_instance: ExecutionEntity = self.child_execution_entities[self.decode_layer_index(relative_locator[0])]
                return located_instance.locate(relative_locator[1:])
            except:
                return None

class IterativeWorkUnit(ControlWorkUnit):
    """
    A base class for managing and executing workflows multiple times or in parallel based on an iterable object.
    
    IterativeWorkUnit serves as an abstract base class that encapsulates common functionalities
    needed for executing a series of workflows iteratively. This class is designed to be extended by
    specific types of work units that require repeated or parallel execution over a set of parameters or data points,
    such as processing parallel tasks in a compute cluster or iterating over a collection of data within a loop.
    
    Features include:
    - Initialization from an iterable object that dictates the number or nature of iterations.
    - Common setup and teardown procedures for each iteration.
    - Abstract methods that must be implemented by subclasses to define specific execution logic.

    This class is intended for use where tasks need to be executed repetitively with variations in input
    or configuration, typical in scenarios like batch processing or parametric studies.

    Attributes:
        sub_workflows (List[WorkFlow]): A list of sub-workflows.
        Other inherited attributes...
        iterable_data (Iterable): The data over which the loop iterates.
        Other inherited attributes...

    Methods:
        __init__: Initializes a IterativeWorkUnit instance.
        from_dict: Creates a IterativeWorkUnit instance from a dictionary configuration.
        reload: Reloads the IterativeWorkUnit instance with updated configuration.
        execute: Executes the loop, iterating over each sub-workflow.
    """
    iterable_data: Iterable
    sub_workflows: List[WorkFlow]

    def __init__(self, unit_dict: dict, workflow: WorkFlow = None, working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), sqlite_filepath: str = str(), debug: bool = False):
        """Initializes an IterativeWorkUnit instance.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (list, optional): A list of the hierarchical indexes to locate the current workunit.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        super().__init__(unit_dict=unit_dict, workflow=workflow, working_directory=working_directory, locator=locator, sqlite_filepath=sqlite_filepath, debug=debug)
        self.sub_workflows = list()
        self.child_execution_entities = self.sub_workflows

    @classmethod
    def from_dict(cls, unit_dict: dict, placeholder_api: Callable, iterative_body_datum_label: str, 
                workflow: WorkFlow = None, working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), sqlite_filepath: str = str(), debug: bool = False) -> IterativeWorkUnit:
        """Initializes an instance of IterativeWorkUnit from a given dictionary.

        As with a normal WorkUnit, the dictionary should contain keys such as `api`, `store_as`, 
        and `args`, which specify the API to be called, the return key to be used for storing the results, 
        and the arguments to the API call, respectively.
        Here, however, the value of `api` must be `loop` to indicate that this WorkUnit is a IterativeWorkUnit, 
        and the arguments contained in `args` are the iteration object and the loop body, where the loop 
        body is treated as a sub-WorkFlow.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            placeholder_api (Callable): a placeholder API marking the necessary arguments for initializing a specific class derived from IterativeWorkUnit.
            iterative_body_datum_label (str): Key indicating the element iterated from ITERABLE_DATA in each sub-workflow.
            workflow (WorkFlow, optional): The parent workflow of this unit.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (list, optional): A list of the hierarchical indexes to locate the current workunit.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            IterativeWorkUnit: An instance of IterativeWorkUnit initialized with the given configuration.
        """
        # Initialize the class.
        unit = cls(unit_dict, workflow, working_directory, locator, sqlite_filepath, debug)
        unit.api = placeholder_api
        unit.status = StatusCode.INITIALIZING

        # Self Inspection.
        unit.self_inspection_and_args_reassembly()
        unit.check_initialization_status()

        # Initialize the Iterative Body.
        data_mapper_to_inherite = unit.generate_data_mapper_to_inherite()

        # Some data should be passed before the initialization of a sub-WorkFlow.
        iterative_body_datum_varname = VarFormatter.unformat_variable(unit.args_dict_to_pass.get(iterative_body_datum_label))
        data_mapper_for_initialize_as_intermediate = {iterative_body_datum_varname: Placeholder.PLACEHOLDER_STR_VALUE}

        # Initialize a PlaceHolder Sub-WorkFlow to perform self-inspection.
        flow_locator = unit.locator.copy()
        flow_locator.append(unit.encode_layer_index())
        placeholder_sub_workflow = WorkFlow.from_list(unit_dict_list=unit.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL),
            debug=unit.debug, intermediate_data_mapper=data_mapper_for_initialize_as_intermediate,
            inherited_data_mapper=data_mapper_to_inherite, locator=flow_locator, control_workunit=unit)
        unit.status = placeholder_sub_workflow.status
        
        # Add error message from placeholder workflow initialization.
        if (placeholder_sub_workflow.error_msg_list):
            unit.error_msg_list += placeholder_sub_workflow.error_msg_list

        return unit
    
    def reload(self, unit_dict: dict, placeholder_api: Callable, iterable_data_label: str, iterative_body_datum_label: str):
        """
        Reloads the IterativeWorkUnit instance with updated configuration.

        This method is used when a user wants to reload the entire task and checks 
        if there are any changes in the parameter values. If the WorkUnit's status 
        indicates an error or pause, it updates the WorkUnit's properties with the 
        new configuration. Otherwise, it compares the current and new configurations 
        to check for updates.

        Args:
            unit_dict (dict): The new configuration dictionary for the WorkUnit.
            placeholder_api (Callable): a placeholder API marking the necessary arguments for initializing a specific class derived from IterativeWorkUnit.
            iterable_data_label (str): Key indicating the object to iterate over in an IterativeWorkUnit.
            iterative_body_datum_label (str): Key indicating the element iterated from ITERABLE_DATA in each sub-workflow.

        Notes:
            - Compares `args_dict_input`, `args_dict_to_pass`, and `params_to_assign_at_execution`
              between the current WorkUnit and a virtual unit created from the new configuration.
            - Also compares `return_key` values.
            - If any of these values are different, sets the status to `StatusCode.SUSPECIOUS_UPDATES`.

        Returns:
            The created virtual WorkUnit instance for further comparison.
        """
        iterative_body_datum_varname_backup = self.args_dict_to_pass.get(iterative_body_datum_label)

        virtual_unit = super().reload(unit_dict=unit_dict, api=placeholder_api)
        reloaded_iterative_data = virtual_unit.args_dict_to_pass.get(iterable_data_label)
        iterative_body_datum_varname = virtual_unit.args_dict_to_pass.get(iterative_body_datum_label)

        data_mapper_to_inherite = self.generate_data_mapper_to_inherite()

        if (iterative_data:=self.args_dict_to_pass_backup.get(iterable_data_label)) == reloaded_iterative_data:
            has_suspecious_updates = False
            has_unresolved_error_or_pause = False

            for iterative_index, iterative_datum in enumerate(iterative_data):
                workflow = self.sub_workflows[iterative_index]
                
                data_mapper_for_update = {iterative_body_datum_varname: iterative_datum}
                data_mapper_for_update.update(data_mapper_to_inherite)

                del workflow.intermediate_data_mapper[iterative_body_datum_varname_backup]
                workflow.intermediate_data_mapper.update(data_mapper_for_update)
                unit_dict_list = virtual_unit.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL)
                workflow.reload(unit_dict_list=unit_dict_list)
                if workflow.status in StatusCode.unexecutable_statuses:
                    self.status = StatusCode.FAILED_INITIALIZATION
                elif workflow.status == StatusCode.SUSPECIOUS_UPDATES:
                    has_suspecious_updates = True
                elif workflow.status in StatusCode.error_or_pause_statuses:
                    has_unresolved_error_or_pause = True
                continue
            
            if self.status not in StatusCode.unexecutable_statuses:
                if has_unresolved_error_or_pause:
                    return
                elif has_suspecious_updates:
                    self.status = StatusCode.SUSPECIOUS_UPDATES
                else:
                    self.status = StatusCode.EXIT_OK
        else:
            _LOGGER.warning(f"Inconsistent {iterable_data_label} in {self.identifier}, all the workflows are to be cleared.")
            self.sub_workflows.clear()
            self.status = StatusCode.SUSPECIOUS_UPDATES

    def update_status_code(self, workflow: WorkFlow) -> None:
        """Update the status code of IterativeWorkUnit (and its derived classes) instance.
        
        Args:
            workflow (WorkFlow): The WorkFlow instance from iterative bodies.
        """
        if self.status == StatusCode.EXIT_WITH_ERROR_AND_PAUSE:
            return
        elif workflow.status == StatusCode.EXIT_WITH_ERROR_AND_PAUSE:
            self.status = StatusCode.EXIT_WITH_ERROR_AND_PAUSE
            return
        elif workflow.status in StatusCode.error_excluding_pause_statuses:
            if self.status in StatusCode.pause_including_error_statuses:
                self.status = StatusCode.EXIT_WITH_ERROR_AND_PAUSE
            else:
                self.status = StatusCode.EXIT_WITH_ERROR_IN_INNER_UNITS
            return
        elif workflow.status in StatusCode.pause_excluding_error_statuses:
            if self.status in StatusCode.error_including_pause_statuses:
                self.status = StatusCode.EXIT_WITH_ERROR_AND_PAUSE
            else:
                self.status = StatusCode.RUNNING_WITH_PAUSE_IN_INNER_UNITS
            return
    
class LoopWorkUnit(IterativeWorkUnit):
    """This class is used to represent a special kind of workunit in a procedural workflow that is used to carry a loop.
    In this kind of workunit, all single units contained in the loop body is packaged together as a subworkflow.
    This class is derived from WorkUnit class.

    - This class encapsulates a loop structure within a workflow, where each iteration of the loop is treated as a sub-workflow.
    - It initializes and executes these sub-workflows successively. Failures in individual sub-workflows do not halt the entire loop's execution.

    Attributes:
        Other inherited attributes...

    Methods:
        loop_unit_placeholder_api: A static placeholder API for initialization.
        __init__: Initializes a LoopWorkUnit instance.
        from_dict: Creates a LoopWorkUnit instance from a dictionary configuration.
        reload: Reloads the LoopWorkUnit instance with updated configuration.
        execute: Executes the loop, iterating over each sub-workflow.
    """

    @staticmethod
    def loop_unit_placeholder_api(workunits: list, loop_data: Iterable, loop_datum_varname: str):
        """Serves as a placeholder API for initializing LoopWorkUnit. This API is not intended for actual execution but marks
        the necessary arguments for initialization.

        Args:
            workunits (dict): A list containing the WorkUnit configuration dicts.
            loop_data (Iterable): The data to be iterated over in the loop. 
                                  This parameter name needs to be consistent with the value of `LOOP_ITERABLE_DATA_LABEL`.
            loop_datum_varname (str): Name of the variable representing each element iterated from ITERABLE_DATA in each workunit of the loop body. 
                                      This parameter name needs to be consistent with the value of `LOOP_BODY_DATUM_LABEL`.
        """
        pass

    def __init__(self, unit_dict: dict, workflow: WorkFlow = None, working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), sqlite_filepath: str = str(), debug: bool = False):
        """Initializes a LoopWorkUnit instance.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (list, optional): A list of the hierarchical indexes to locate the current workunit.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        super().__init__(unit_dict=unit_dict, workflow=workflow, working_directory=working_directory, locator=locator, sqlite_filepath=sqlite_filepath, debug=debug)
        return

    @classmethod
    def from_dict(cls, unit_dict: dict, workflow: WorkFlow = None, working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), sqlite_filepath: str = str(), debug: bool = False) -> LoopWorkUnit:
        """Initializes an instance of LoopWorkUnit from a given dictionary.

        As with a normal WorkUnit, the dictionary should contain keys such as `api`, `store_as`, 
        and `args`, which specify the API to be called, the return key to be used for storing the results, 
        and the arguments to the API call, respectively.
        Here, however, the value of `api` must be `loop` to indicate that this WorkUnit is a LoopWorkUnit, 
        and the arguments contained in `args` are the iteration object and the loop body, where the loop 
        body is treated as a sub-WorkFlow.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            workflow (WorkFlow, optional): The parent workflow of this unit.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (list, optional): A list of the hierarchical indexes to locate the current workunit.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            LoopWorkUnit: An instance of LoopWorkUnit initialized with the given configuration.
        """
        _LOGGER.info(f"Initializing a {__class__.__name__} now...")
        unit = super().from_dict(unit_dict=unit_dict, placeholder_api=__class__.loop_unit_placeholder_api, 
                            iterative_body_datum_label=LOOP_BODY_DATUM_LABEL, workflow=workflow, working_directory=working_directory,
                            locator=locator, sqlite_filepath=sqlite_filepath, debug=debug)
        return unit

    def reload(self, unit_dict: dict):
        """
        Reloads the LoopWorkUnit instance with updated configuration.

        This method is used when a user wants to reload the entire task and checks 
        if there are any changes in the parameter values. If the WorkUnit's status 
        indicates an error or pause, it updates the WorkUnit's properties with the 
        new configuration. Otherwise, it compares the current and new configurations 
        to check for updates.

        Args:
            unit_dict (dict): The new configuration dictionary for the WorkUnit.

        Notes:
            - Compares `args_dict_input`, `args_dict_to_pass`, and `params_to_assign_at_execution`
              between the current WorkUnit and a virtual unit created from the new configuration.
            - Also compares `return_key` values.
            - If any of these values are different, sets the status to `StatusCode.SUSPECIOUS_UPDATES`.

        Returns:
            The created virtual WorkUnit instance for further comparison.
        """
        super().reload(unit_dict=unit_dict, placeholder_api=__class__.loop_unit_placeholder_api, 
                    iterable_data_label=LOOP_ITERABLE_DATA_LABEL, iterative_body_datum_label=LOOP_BODY_DATUM_LABEL)

    def execute(self) -> Tuple[str, Any]:
        """Executes the loop, iterating over each sub-workflow in `sub_workflows`.

        In each iteration, the corresponding sub-workflow is executed, and its results are collected.
        Failures in sub-workflows do not stop the execution of the entire loop.

        Returns:
            tuple: A tuple containing the return key and a list of results from each loop iteration.
        """
        if self.status in StatusCode.skippable_statuses:
            return self.return_key, self.return_value

        self.return_value: dict = dict()
        data_mapper_to_inherite = self.generate_data_mapper_to_inherite()
        self.update_args_value_from_data_mapper()
        _LOGGER.info(f'We are about to execute {self.__class__.__name__} {self.identifier} ...')
        self.status = StatusCode.RUNNING
        loop_body_datum_varname = self.args_dict_to_pass.get(LOOP_BODY_DATUM_LABEL)
        loop_data = self.args_dict_to_pass.get(LOOP_ITERABLE_DATA_LABEL)

        is_existing_subflows = (len(loop_data) == len(self.sub_workflows))
        if (not is_existing_subflows):  # If we don't use existing sub-workflows, clear all the sub-workflows.
            self.sub_workflows.clear()

        # Add sub_workflow(s) to the execution queue.
        for loop_index, loop_datum in enumerate(loop_data):
            workflow: WorkFlow
            data_mapper_for_initialization = {loop_body_datum_varname: loop_datum}
            
            if (is_existing_subflows):
                workflow = self.sub_workflows[loop_index]
                workflow.intermediate_data_mapper.update(data_mapper_for_initialization)
                workflow.inherited_data_mapper.update(data_mapper_to_inherite)
            else:
                flow_locator = self.locator.copy()
                flow_locator.append(self.encode_layer_index(loop_index))
                workflow = WorkFlow.from_list(
                    unit_dict_list=self.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL), 
                    working_directory=self.working_directory, debug=self.debug, 
                    intermediate_data_mapper=data_mapper_for_initialization,
                    inherited_data_mapper=data_mapper_to_inherite, locator=flow_locator,
                    control_workunit=self, sqlite_filepath=self.sqlite_filepath)
                self.sub_workflows.append(workflow)
            continue

        # Execute the sub-workflows one by one.
        for loop_index, workflow in enumerate(self.sub_workflows):
            _LOGGER.debug(f'Executing {workflow.identifier} ...')

            # Temporarily change directory before execution.
            flow_working_subfolder = f"{VarFormatter.unformat_variable(self.args_dict_to_pass.get(LOOP_BODY_DATUM_LABEL))}_{loop_index}"
            safe_mkdir(flow_working_subfolder)
            chdir(flow_working_subfolder)

            workflow_return = workflow.execute()

            # Change directory back after execution.
            chdir("..")

            self.return_value[self.encode_layer_index(loop_index)] = workflow_return
            self.update_status_code(workflow=workflow)
            continue

        if (self.status == StatusCode.RUNNING):
            self.status = StatusCode.EXIT_OK
        return self.return_key, self.return_value

class ClusterBatchWorkUnit(IterativeWorkUnit):
    """Manages and orchestrates computational tasks across a cluster computing environment.
    
    The ClusterBatchWorkUnit class is derived from ControlWorkUnit and is designed to handle the 
    distribution and execution of computational tasks on multiple nodes within a computing cluster.
    This class facilitates the submission, tracking, and scheduling of parallel tasks, leveraging
    the power of cluster resources to perform large-scale computations more efficiently.
    
    - This class encapsulates a cluster job within a workflow, where each instance of the workflow is treated as a sub-workflow.
    - It initializes and submits these sub-workflows to computational cluster. Failures in individual sub-workflows do not halt the execution of the entire WorkUnit.

    Attributes:
        iterable_data (Iterable): The data over which the WorkUnit encapsulates cluster batch jobs.
        max_simultaeneous_jobs (int): The maximum number of simultaneous tasks submitted to the cluster.
        Other inherited attributes...

    Methods:
        cluster_batch_unit_placeholder_api: A static placeholder API for initialization.
        __init__: Initializes a ClusterBatchWorkUnit instance.
        from_dict: Creates a ClusterBatchWorkUnit instance from a dictionary configuration.
        execute: Executes the cluster batch, iterating over each sub-workflow.
    """
    iterable_data: Iterable
    max_simultaeneous_jobs: int

    @staticmethod
    def cluster_batch_placeholder_api(workunits: list, batch_data: Iterable, batch_datum_varname: str, max_simultaeneous_jobs: int = DEFAULT_CLUSTER_JOB_CAPABILITY):
        """Serves as a placeholder API for initializing ClusterBatchWorkUnit. This API is not intended for actual execution but marks
        the necessary arguments for initialization.

        Args:
            workunits (dict): A list containing the WorkUnit configuration dicts.
            batch_data (Iterable): The data to be iterated over in the loop. 
                                  This parameter name needs to be consistent with the value of `BATCH_ITERABLE_DATA_LABEL`.
            batch_datum_varname (str): Name of the variable representing each element iterated from ITERABLE_DATA in each workunit of the loop body. 
                                      This parameter name needs to be consistent with the value of `BATCH_BODY_DATUM_LABEL`.
            max_simultaeneous_jobs (int, optional): Due to the resource limit of the Cluster, we have to set a maximum number of simultaneous tasks 
                                    submitted to the cluster (Default: the value of DEFAULT_CLUSTER_JOB_CAPABILITY).
        """
        pass

    def __init__(self, unit_dict: dict, workflow: WorkFlow = None, working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), sqlite_filepath: str = str(), debug: bool = False):
        """Initializes a ClusterBatchWorkUnit instance.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (list, optional): A list of the hierarchical indexes to locate the current workunit.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        super().__init__(unit_dict=unit_dict, workflow=workflow, working_directory=working_directory, locator=locator, sqlite_filepath=sqlite_filepath, debug=debug)
        return
    
    @classmethod
    def from_dict(cls, unit_dict: dict, workflow: WorkFlow = None, working_directory: str = BASE_DIRECTORY, locator: List[int] = list(), sqlite_filepath: str = str(), debug: bool = False) -> ClusterBatchWorkUnit:
        """Initializes an instance of ClusterBatchWorkUnit from a given dictionary.

        As with a normal WorkUnit, the dictionary should contain keys such as `api`, `store_as`, 
        and `args`, which specify the API to be called, the return key to be used for storing the results, 
        and the arguments to the API call, respectively.
        Here, however, the value of `api` must be `cluster_batch` to indicate that this WorkUnit is a ClusterBatchWorkUnit, 
        and the arguments contained in `args` are the iteration object and the batch body, where the loop 
        body is treated as a sub-WorkFlow.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            workflow (WorkFlow, optional): The parent workflow of this unit.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            locator (list, optional): A list of the hierarchical indexes to locate the current workunit.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            ClusterBatchWorkUnit: An instance of ClusterBatchWorkUnit initialized with the given configuration.
        """
        _LOGGER.info(f"Initializing a {__class__.__name__} now...")
        unit: __class__ = super().from_dict(unit_dict=unit_dict, placeholder_api=__class__.cluster_batch_placeholder_api, 
                            iterative_body_datum_label=BATCH_BODY_DATUM_LABEL, workflow=workflow, working_directory=working_directory, locator=locator, 
                            sqlite_filepath=sqlite_filepath, debug=debug)
        unit.max_simultaeneous_jobs = unit.args_dict_to_pass.get("max_simultaeneous_jobs", DEFAULT_CLUSTER_JOB_CAPABILITY)
        return unit

    def reload(self, unit_dict: dict):
        """
        Reloads the ClusterBatchWorkUnit instance with updated configuration.

        This method is used when a user wants to reload the entire task and checks 
        if there are any changes in the parameter values. If the WorkUnit's status 
        indicates an error or pause, it updates the WorkUnit's properties with the 
        new configuration. Otherwise, it compares the current and new configurations 
        to check for updates.

        Args:
            unit_dict (dict): The new configuration dictionary for the WorkUnit.

        Notes:
            - Compares `args_dict_input`, `args_dict_to_pass`, and `params_to_assign_at_execution`
              between the current WorkUnit and a virtual unit created from the new configuration.
            - Also compares `return_key` values.
            - If any of these values are different, sets the status to `StatusCode.SUSPECIOUS_UPDATES`.

        Returns:
            The created virtual WorkUnit instance for further comparison.
        """
        super().reload(unit_dict=unit_dict, placeholder_api=__class__.cluster_batch_placeholder_api, 
                    iterable_data_label=BATCH_ITERABLE_DATA_LABEL, iterative_body_datum_label=BATCH_BODY_DATUM_LABEL)

    def execute(self) -> Tuple[str, Any]:
        """Executes the cluster batch, iterating over each sub-workflow in `sub_workflows`.

        In each iteration, the corresponding sub-workflow is executed, and its results are collected.
        Failures in sub-workflows do not stop the execution of the entire loop.

        Returns:
            tuple: A tuple containing the return key and a list of results from each loop iteration.
        """
        if self.status in StatusCode.skippable_statuses:
            return self.return_key, self.return_value

        self.return_value: dict = dict()
        data_mapper_to_inherite = self.generate_data_mapper_to_inherite()
        self.update_args_value_from_data_mapper()
        _LOGGER.info(f'We are about to execute {self.__class__.__name__} {self.identifier} ...')
        self.status = StatusCode.RUNNING
        batch_body_datum_varname = self.args_dict_to_pass.get(BATCH_BODY_DATUM_LABEL)
        batch_data = self.args_dict_to_pass.get(BATCH_ITERABLE_DATA_LABEL)

        is_existing_subflows = (len(batch_data) == len(self.sub_workflows))
        if (not is_existing_subflows):  # If we don't use existing sub-workflows, clear all the sub-workflows.
            self.sub_workflows.clear()

        # Add sub_workflow(s) to the execution queue.
        for batch_index, batch_datum in enumerate(batch_data):
            workflow: WorkFlow
            data_mapper_for_initialization = {batch_body_datum_varname: batch_datum}
            
            if (is_existing_subflows):
                workflow = self.sub_workflows[batch_index]
                workflow.intermediate_data_mapper.update(data_mapper_for_initialization)
                workflow.inherited_data_mapper.update(data_mapper_to_inherite)
            else:
                flow_locator = self.locator.copy()
                flow_locator.append(self.encode_layer_index(batch_index))
                workflow = WorkFlow.from_list(
                    unit_dict_list=self.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL),
                    working_directory=self.working_directory, debug=self.debug, 
                    intermediate_data_mapper=data_mapper_for_initialization,
                    inherited_data_mapper=data_mapper_to_inherite, locator=flow_locator,
                    control_workunit=self, sqlite_filepath=self.sqlite_filepath)
                self.sub_workflows.append(workflow)
            continue

        # Execute the sub-workflows one by one. TODO: 这里还是从 Loop 复制过来的状态，需要重新写，涉及 sbatch 提交等事务。
        for batch_index, workflow in enumerate(self.sub_workflows):
            _LOGGER.debug(f'Executing {workflow.identifier} ...')

            # Temporarily change directory before execution.
            current_dir = getcwd()
            flow_working_subfolder = f"{VarFormatter.unformat_variable(self.args_dict_to_pass.get(LOOP_BODY_DATUM_LABEL))}_{batch_index}"
            safe_mkdir(flow_working_subfolder)
            chdir(flow_working_subfolder)

            workflow_return = workflow.execute()

            # Change directory back after execution.
            chdir(current_dir)

            self.return_value[self.encode_layer_index(batch_index)] = workflow_return
            self.update_status_code(workflow=workflow)
            continue

        if (self.status == StatusCode.RUNNING):
            self.status = StatusCode.EXIT_OK
        return self.return_key, self.return_value

#endregion

#region Entry Point.

class GeneralWorkUnit(ControlWorkUnit):
    """This class is used to represent a special kind of unit of work that would normally be used at the outermost level 
    of a workflow configuration, existing as a carrier for the entire workflow.
    This class is derived from WorkUnit class.

    Attributes:
        sub_workflow (WorkFlow): The workflow it carries.
        working_directory (str): Indicates where to run the job. Default to current working directory.
        save_snapshot (bool): Whether to automatically save the status as a pickle file 
                            when it exits due to Error, Pause, or Completion. Default False.
        data_mapper_for_init (dict): An data mapper for initialization, which can be used to pass in data of types that are inconvenient to record in JSON.
        Other inherited attributes...

    Methods:
        general_unit_placeholder_api: A static placeholder API for initialization.
        __init__: Initializes a GeneralWorkUnit instance.
        from_dict: Creates a GeneralWorkUnit instance from a dictionary configuration.
        execute: Executes the GeneralWorkUnit, i.e., executes the sub_workflow.
    """
    sub_workflow: WorkFlow
    sqlite_filename: str
    save_snapshot: bool
    latest_pickle_filepath: str
    data_mapper_for_init: dict

    @staticmethod
    def general_unit_placeholder_api(workunits: list, **kwargs):
        """Serves as a placeholder API for initializing LoopWorkUnit. This API is not intended for actual execution but marks
        the necessary arguments for initialization.

        Args:
            workunits (list): A list containing the WorkUnit configuration dicts.
        """
        return

    def __init__(self, unit_dict: dict, 
                workflow: WorkFlow = None, 
                working_directory: str = BASE_DIRECTORY, 
                save_snapshot: bool = False, 
                sqlite_filepath: str = str(), 
                debug: bool = False):
        """Initializes a GeneralWorkUnit instance.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to. Usually None.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            save_snapshot (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            sqlite_filepath (str): The filepath to your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        super().__init__(unit_dict=unit_dict, workflow=workflow, working_directory=working_directory, sqlite_filepath=sqlite_filepath, debug=debug)
        self.save_snapshot = save_snapshot
        self.api_key = GENERAL_API_KEY
        self.latest_pickle_filepath = str()
        self.status = StatusCode.CREATED

    @classmethod
    def from_dict(cls, unit_dict: dict, working_directory: str = BASE_DIRECTORY, 
                sqlite_filename: str = DEFAULT_SQLITE_FILENAME, overwrite_database: bool = False,
                save_snapshot: bool = False, debug: bool = False, data_mapper_for_init: dict = dict()) -> GeneralWorkUnit:
        """Initializes an instance of GeneralWorkUnit from a given dictionary.

        Here, however, `api` must be `general` to indicate that the WorkUnit is a GeneralWorkUnit, 
        and the arguments contained in `args` are an assemble of WorkUnits forming a WorkFlow and 
        required information to execute the GeneralWorkUnit (the whole workflow).

        Args:
            unit_dict (dict): Configuration dictionary parsed from Json configuration containing information to initialize a task.
            working_directory (str, optional): Indicates where to run the job. Default to current directory.
            save_snapshot (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            sqlite_filename (str): The filename of your sqlite database file for persistent storage.
            debug (bool, optional): Indicates whether to run in debug mode.
            data_mapper_for_init (dict, optional): An data mapper for initialization, which can be used to pass in data of types that are inconvenient to record in JSON.

        Returns:
            LoopWorkUnit: An instance of LoopWorkUnit initialized with the given configuration.
        """
        # Handle working directory.
        # Highest priority: use the working directory passed in by the initialization function; 
        # Lowest priority: use BASE_DIRECTORY.
        if path.abspath(working_directory) == path.abspath(BASE_DIRECTORY):
            if (work_dir_from_unit_dict:=unit_dict.get(WORKUNIT_ARGUMENT_LIST_KEY).get(WORKING_DIRECTORY_KEY)):
                working_directory = work_dir_from_unit_dict

        unit_dict[WORKUNIT_ARGUMENT_LIST_KEY][WORKING_DIRECTORY_KEY] = working_directory    # Update working directory to the unit dict.
        safe_mkdir(working_directory)
        sqlite_filepath = path.join(working_directory, sqlite_filename)

        # Temporarily change directory before initialization.
        chdir(working_directory)

        # Initialize the class.
        unit = cls(unit_dict=unit_dict, working_directory=working_directory, save_snapshot=save_snapshot, sqlite_filepath=sqlite_filepath, debug=debug)
        database_session = unit.create_database_session(overwrite_database=overwrite_database)
        database_session.close()
        if not save_snapshot:
            not_save_snapshot_warning_msg = "The `save_snapshot` value is not set to `True`, so this GeneralWorkUnit instance will not automatically save its state on error or pause and continue computing later after corrections."
            _LOGGER.warning(not_save_snapshot_warning_msg)
            if not debug:
                sleep_seconds = 5.0
                _LOGGER.warning(f"Press CTRL+C to cancel the job now, or we continue in {sleep_seconds} seconds...\n")
                sleep(sleep_seconds)

        unit.api = GeneralWorkUnit.general_unit_placeholder_api
        unit.data_mapper_for_init = data_mapper_for_init

        # Self Inspection.
        unit.self_inspection_and_args_reassembly()
        unit.check_initialization_status()
        if (unit.status == StatusCode.READY_TO_START):
            unit.status = StatusCode.INITIALIZING
        
        # Initialize the WorkFlow.
        data_mapper_for_inherit = unit.args_dict_to_pass.copy()
        data_mapper_for_inherit.update(unit.data_mapper_for_init)
        
        del data_mapper_for_inherit[CONTROL_BODY_WORKUNITS_LABEL]  # Pass kwargs to the workflow.

        if (unit_dict_list:=unit.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL)):
            unit.sub_workflow = WorkFlow.from_list(unit_dict_list, debug=debug, inherited_data_mapper=data_mapper_for_inherit, control_workunit=unit, sqlite_filepath=unit.sqlite_filepath)
            if unit.sub_workflow.status == StatusCode.READY_TO_START:
                unit.status = StatusCode.READY_TO_START
            else:
                unit.status = StatusCode.FAILED_INITIALIZATION
        else:
            error_msg = "Initializing GeneralWorkUnit with empty `workunits`. WorkUnit(s) are expected."
            unit.status = StatusCode.FAILED_INITIALIZATION
            _LOGGER.error(error_msg)
            raise ValueError(error_msg)
    
        # Change directory back after initialization.
        chdir(BASE_DIRECTORY)
        return unit
    
    @classmethod
    def from_json_string(cls, json_str: str, 
                        working_directory: str = BASE_DIRECTORY, save_snapshot: bool = False, 
                        sqlite_filename: str = DEFAULT_SQLITE_FILENAME, overwrite_database: bool = False,
                        debug: bool = False, data_mapper_for_init: dict = dict()) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a serialized json string.
        
        Args:
            json_str (str): The serialized json string.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            save_snapshot (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            sqlite_filename (str, optional): The filename of the sqlite database file in your working directory.
            overwrite_database: If overwrite the database file when it exists. Default False.
            debug (bool, optional): Indicates whether to run in debug mode. Default False.
            data_mapper_for_init (dict, optional): An data mapper for initialization, which can be used to pass in data of types that are inconvenient to record in JSON.
        
        Returns:
            An instance of WorkFlow.
        """
        unit_dict = loads(json_str)
        return GeneralWorkUnit.from_dict(unit_dict=unit_dict, working_directory=working_directory, 
                                        save_snapshot=save_snapshot, sqlite_filename=sqlite_filename, 
                                        overwrite_database=overwrite_database, debug=debug, 
                                        data_mapper_for_init=data_mapper_for_init)

    @classmethod
    def from_json_file_object(cls, json_fobj: TextIOWrapper, 
                            working_directory: str = BASE_DIRECTORY, save_snapshot: bool = False, 
                            sqlite_filename: str = DEFAULT_SQLITE_FILENAME, overwrite_database: bool = False,
                            debug: bool = False, data_mapper_for_init: dict = dict()) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a json file object.
        
        Args:
            json_fobj (TextIOWrapper): The instance of a json file containing GeneralWorkUnit configuration information.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            save_snapshot (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            sqlite_filename (str, optional): The filename of the sqlite database file in your working directory.
            overwrite_database: If overwrite the database file when it exists. Default False.
            debug (bool, optional): Indicates whether to run in debug mode. Default False.
            data_mapper_for_init (dict, optional): An data mapper for initialization, which can be used to pass in data of types that are inconvenient to record in JSON.
        
        Returns:
            An instance of GeneralWorkUnit.
        """
        json_str = json_fobj.read()
        return GeneralWorkUnit.from_json_string(json_str=json_str, working_directory=working_directory, 
                                                save_snapshot=save_snapshot, sqlite_filename=sqlite_filename, overwrite_database=overwrite_database,
                                                debug=debug, data_mapper_for_init=data_mapper_for_init)

    @classmethod
    def from_json_filepath(cls, json_filepath: str, 
                        working_directory: str = BASE_DIRECTORY, save_snapshot: bool = False, 
                        sqlite_filename: str = DEFAULT_SQLITE_FILENAME, overwrite_database: bool = False,
                        debug: bool = False, data_mapper_for_init: dict = dict()) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a json filepath.
        
        Args:
            json_filepath (str): The path to a json file containing GeneralWorkUnit configuration information.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            save_snapshot (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            sqlite_filename (str, optional): The filename of the sqlite database file in your working directory.
            overwrite_database: If overwrite the database file when it exists. Default False.
            debug (bool, optional): Indicates whether to run in debug mode. Default False.
            data_mapper_for_init (dict, optional): An data mapper for initialization, which can be used to pass in data of types that are inconvenient to record in JSON.
        
        Returns:
            An instance of WorkFlow.
        """
        with open(json_filepath) as fobj:
            return GeneralWorkUnit.from_json_file_object(fobj, working_directory=working_directory, 
                                                    save_snapshot=save_snapshot, sqlite_filename=sqlite_filename, 
                                                    overwrite_database=overwrite_database, debug=debug, 
                                                    data_mapper_for_init=data_mapper_for_init)
        
    def execute(self) -> Tuple[str, Any]:
        """Executes the GeneralWorkUnit.

        Returns:
            tuple: A tuple containing the return key and the intermediate data from the workflow.
        """
        if not self.is_executable:
            _LOGGER.error(FAILED_INITIALIZATION_ERROR_MSG)
            raise ValueError(self)
        else:
            _LOGGER.info(SUCCESS_INITIALIZATION_MSG)
            self.status = StatusCode.RUNNING

            # Temporarily change directory before execution.
            chdir(self.working_directory)
            self.return_value = self.sub_workflow.execute()

            # Change directory back after execution.
            chdir(BASE_DIRECTORY)

            if (self.sub_workflow.status in StatusCode.error_excluding_pause_statuses):
                self.status = StatusCode.EXIT_WITH_ERROR_IN_INNER_UNITS
            elif ((self.sub_workflow.status in StatusCode.pause_excluding_error_statuses) and (self.status not in StatusCode.error_excluding_pause_statuses)):
                self.status = StatusCode.RUNNING_WITH_PAUSE_IN_INNER_UNITS
            else:
                self.status = StatusCode.EXIT_OK
            if (self.save_snapshot):
                pickle_filepath = self.dump_snapshot_file()
                self.latest_pickle_filepath = pickle_filepath
                _LOGGER.info(f"The current running status has been saved to '{pickle_filepath}'. Please go to this location to view your file(s).")
            return self.return_key, self.return_value
   
    def locate(self, locator: list) -> ExecutionEntity:
        """Recursive locating of ExecutionEntity instance based on their locators.

        Note: This method will return None for an unexecuted ExecutionEntity instance that is inside a ControlWorkUnit instance other than GeneralWorkUnit.
        
        Args:
            locator (list): The locator of the target ExecutionEntity instance.

        Returns:
            The located ExecutionEntity instance. None if nothing is located.
        """
        if len(locator) == 0:
            return self
        else:
            try:
                return self.sub_workflow.locate(locator)
            except:
                return None
    
    def reload(self, unit_dict: dict, data_mapper_for_reload: dict = dict()) -> None:
        """Reload task from a GeneralWorkUnit dict.
        
        Note: If the architecture is changed, the task is not able to be reloaded.

        Args:
            unit_dict (dict): Modified configuration dictionary parsed from Json configuration containing information to initialize a task.
            data_mapper_for_reload (dict): An data mapper for updating the initialization data mapper (data_mapper_for_init).
        """
        architecture = __class__.read_architecture(unit_dict=unit_dict)
        if (architecture != self.architecture):
            raise ValueError("Unable to reload the GeneralWorkUnit. We currently do not support reloading the task after architectural changes.")

        # Temporarily change directory before initialization.
        chdir(self.working_directory)
        
        virtual_unit = super().reload(unit_dict=unit_dict, api=__class__.general_unit_placeholder_api)

        if (not set(data_mapper_for_reload.items()).issubset(set(self.data_mapper_for_init.items()))):
            # If not all the items in `data_mapper_for_reload` is already in `self.data_mapper_for_init`, then updates exist.
            self.data_mapper_for_init.update(data_mapper_for_reload)
            self.status = StatusCode.SUSPECIOUS_UPDATES if self.status != StatusCode.FAILED_INITIALIZATION else StatusCode.FAILED_INITIALIZATION
        
        if (self.status == StatusCode.FAILED_INITIALIZATION):
            raise ValueError("Unable to reload the GeneralWorkUnit. Format errors exist your json configuration file!")
        elif (self.status == StatusCode.SUSPECIOUS_UPDATES):
            self.data_mapper_for_init.update(data_mapper_for_reload)

            # Reload the workflow.
            unit_dict_list = virtual_unit.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL)

            # Pass data mapper to the workflow.
            data_mapper_for_inherit = virtual_unit.args_dict_to_pass.copy()
            data_mapper_for_inherit.update(self.data_mapper_for_init)
            del data_mapper_for_inherit[CONTROL_BODY_WORKUNITS_LABEL]

            self.sub_workflow.inherited_data_mapper.update(data_mapper_for_inherit)
            self.sub_workflow.reload(unit_dict_list=unit_dict_list)

        # Change directory back after initialization.
        chdir(BASE_DIRECTORY)

    def reload_json_string(self, json_str: str, data_mapper_for_reload: dict = dict()) -> None:
        """Reload task from a serialized json string.
        
        Args:
            json_str (str): The serialized json string.
            data_mapper_for_reload (dict): An data mapper for updating the initialization data mapper (data_mapper_for_init).
        """
        unit_dict = loads(json_str)
        self.reload(unit_dict, data_mapper_for_reload=data_mapper_for_reload)
        return

    def reload_json_file_object(self, json_fobj: TextIOWrapper, data_mapper_for_reload: dict = dict()) -> None:
        """Reload task from JSON format configuration file object.
        
        Args:
            json_fobj (str): The modified json file object.
            data_mapper_for_reload (dict): An data mapper for updating the initialization data mapper (data_mapper_for_init).
        """
        json_str = json_fobj.read()
        self.reload_json_string(json_str, data_mapper_for_reload=data_mapper_for_reload)
        return

    def reload_json_filepath(self, json_filepath: str, data_mapper_for_reload: dict = dict()) -> None:
        """Reload task from JSON format configuration filepath.
        
        Args:
            json_filepath (str): The modified json filepath.
            data_mapper_for_reload (dict): An data mapper for updating the initialization data mapper (data_mapper_for_init).
        """
        with open(json_filepath) as fobj:
            self.reload_json_file_object(fobj, data_mapper_for_reload=data_mapper_for_reload)
        return

#endregion
