#! python3
# -*- encoding: utf-8 -*-
"""
Workflow Interpreter.

TODO (Zhong): Typing type annotation inspection. (Completed 2024-02-16 23:16 UTC-6)
TODO (Zhong): Annotation inspection in intermediate_data_mapper. (Completed 2024-03-01 01:00 UTC-6)
TODO (Zhong): APIs entitled `loop`. (Completed 2024-02-26 20:01 UTC-6)
TODO (Zhong): How to handle intermediate data transfer between parent and child workflows. Child workflow should be able to read the data in parent flow. (Completed 2024-02-28 06:06 UTC-6)
TODO (Zhong): Continue computing.
TODO (Zhong): How to indicate error from inner units. (Completed 2024-03-13 19:00 UTC-5)
TODO (Zhong): Data Output.

@File    :   workflow.py
@Created :   2024/02/05 22:53
@Author  :   Zhong, Yinjie
@Contact :   yinjie.zhong@vanderbilt.edu
"""

# Here put the import lib.
from __future__ import annotations
from os import getcwd, path
from io import TextIOWrapper, FileIO
from time import sleep
from inspect import _empty, Parameter, signature
from json import load, loads
from typing import Any, Dict, List, Type, Union, Tuple
from collections.abc import Iterable, Callable
from datetime import datetime
import pickle

from enzy_htp.core.logger import _LOGGER

from .config import (
    SCIENCE_API_MAPPER, 
    JsonConfigVariableFormatter as VarFormatter, 
)

DEFAULT_WORK_DIR = getcwd() # Default Working Directory: Current working directory.

# In the Json file used to define the workflow, 
# keys in the mapper are defined by the developer,
# which are not user-definable.
WORKUNIT_API_NAME_KEY = "api"               # Key indicating the API name.
WORKUNIT_RETURN_VALUE_KEY = "store_as"      # Key indicating the key where the return value of the workunit is stored.
WORKUNIT_ARGUMENT_LIST_KEY = "args"         # Key indicates the arguments to be passed into the API.
CONTROL_BODY_WORKUNITS_LABEL = "workunits"  # Key indicating the list of workunits contained in the body in a unit playing a control role.
LOOP_API_KEY = "loop"                       # The key of a loop workunit.
LOOP_ITERABLE_DATA_LABEL = "loop_data"      # Key indicating the object to iterate over in a LoopWorkUnit.
LOOP_BODY_DATUM_LABEL = "loop_datum_varname"# Key indicating the element iterated from ITERABLE_DATA in each workunit of the loop body.
GENERAL_API_KEY = "general"                 # The key of a general workunit.

# API Keys for flow control, such as `loop`, `general`.
CONTROL_API_KEYS = [
    LOOP_API_KEY,
    GENERAL_API_KEY,
]

# Some string values used globally are defined here.
FAILED_INITIALIZATION_ERROR_MSG = "The initialized work has not yet passed the self-inspection, so it is not allowed to be executed!"
SUCCESS_INITIALIZATION_MSG = "The initialized work has successfully passed the self-inspection. Now proceeding to execution."
PLACEHOLDER_VALUE = "The world is yours and ours, but in the final analysis it's you young people's."  # Value for placeholder.

class StatusCode():
    """
    Class representing various execution statuses for workunits and workflows.

    Notes:
        When a WorkFlow instance or WorkUnit instance is marked as `RUNNING_WITH_PAUSE_IN_INNER_UNITS` or `EXPECTED_PAUSE`, 
        it freezes the execution of all subsequent instances of the WorkUnit for the WorkFlow instance in which it resides
        (a loop ControlWorkUnit or parallel ControlWorkUnit freezes execution of the current loop/parallel instance and 
        proceeds directly to the next loop/parallel instance until all instances of the ControlWorkUnit in which it resides 
        have been executed), and then saves the current GeneralWorkUnit instance state to a snapshot pickle file.

    Notes in Chinese Language for Reference (consistent with the meaning of the English version) / 汉语备注，备忘用，与英语文本含义一致:
        当一个 WorkFlow 实例或 WorkUnit 实例被标记为 `RUNNING_WITH_PAUSE_IN_INNER_UNITS` 或 `EXPECTED_PAUSE` 时，
        将冻结其所在的工作流实例的后续所有 WorkUnit 实例的执行（循环式 ControlWorkUnit 或并行式 ControlWorkUnit 
        会冻结当前循环体/并行体实例，直接继续执行下一个循环体/并行体实例，直至其所在的 ControlWorkUnit 中的所有循环体/并行体实例被执行完毕），
        然后保存当前 GeneralWorkUnit 实例状态的快照至 pickle 文件。

    Attributes:
        CREATED (int): The initial status when a WorkUnit or WorkFlow instance is created. Value: -9
        INITIALIZING (int): The status when a WorkUnit or WorkFlow instance is undergoin initialization. Value: -8
        READY_TO_START (int): The status when a WorkUnit or WorkFlow instance has passed the self-inspection but hasn't yet been started. Value: -4
        RUNNING (int): Indicates that the workunit or workflow is currently in execution. Value: -3
        PAUSE_IN_INNER_UNITS (int): Specific to WorkFlow and ControlWorkUnit instances. Indicates
                                    that the workflow is running but an inner unit is paused as expected. Value: -2
        EXPECTED_PAUSE (int): Specific to Basic WorkUnit instances. Indicates that a unit is paused and its outer
                              layers should be marked as `RUNNING_WITH_PAUSE_IN_INNER_UNITS`. Value: -1
        EXIT_WITH_OK (int): Indicates successful completion of the work unit or workflow. Value: 0
        ERROR_IN_INNER_UNITS (int): Specific to WorkFlow and ControlWorkUnit. Indicates error(s) in the
                                    inner units of a workflow. Value: 1
        EXIT_WITH_ERROR (int): Specific to Basic WorkUnit instances. Indicates that the work unit or workflow
                               exited with an error. Value: 2
        EXIT_WITH_ERROR_AND_PAUSE (int): Specific to WorkFlow and ControlWorkUnit. Indicates the coexistence of error(s)
                                        and expected pause(s) in the inner units of a workflow. Value: 3
        FAILED_INITIALIZATION (int): Indicates that the initialization of the workunit or workflow failed. Value: 9
    """
    CREATED = -9
    INITIALIZING = -8
    READY_TO_START = -4
    RUNNING = -3
    RUNNING_WITH_PAUSE_IN_INNER_UNITS = -2   # For WorkFlow and ControlWorkUnit only.
    EXPECTED_PAUSE = -1                 # For Basic WorkUnit only. If a unit is paused, set its outer layers as `RUNNING_WITH_PAUSE_IN_INNER_UNITS`.
    EXIT_OK = 0
    EXIT_WITH_ERROR_IN_INNER_UNITS = 1  # For WorkFlow and ControlWorkUnit only.
    EXIT_WITH_ERROR = 2                 # For Science API only.
    EXIT_WITH_ERROR_AND_PAUSE = 3       # For WorkFlow and ControlWorkUnit only.
    FAILED_INITIALIZATION = 9

    pause_excluding_error_statuses = [EXPECTED_PAUSE, RUNNING_WITH_PAUSE_IN_INNER_UNITS]    # For logical judgment only.
    pause_including_error_statuses = [EXPECTED_PAUSE, RUNNING_WITH_PAUSE_IN_INNER_UNITS, EXIT_WITH_ERROR_AND_PAUSE] # For logical judgment only.
    error_excluding_pause_statuses = [EXIT_WITH_ERROR, EXIT_WITH_ERROR_IN_INNER_UNITS]      # For logical judgment only.
    error_including_pause_statuses = [EXIT_WITH_ERROR, EXIT_WITH_ERROR_IN_INNER_UNITS, EXIT_WITH_ERROR_AND_PAUSE]   # For logical judgment only.
    unexecutable_statuses = [CREATED, INITIALIZING, FAILED_INITIALIZATION]                  # For logical judgment only.

class WorkFlow():
    """Rrepresents a procedural WorkFlow primarily consists of a number of WorkUnit instances, 
    including data mappers, self-inspection flags, to ensure its operation.

    Each WorkFlow Instance will perform self-inspection during initialization. If the self-inspection fails, it cannot be executed.
    
    If an error occurs during the execution of a single WorkUnit instance, the workflow will record the error, 
    jump out of the WorkUnit instance, and execute the next one until the execution of the entire workflow is completed.
    
    Attributes:
        debug (bool): Flag indicating whether debugging mode is enabled.
        error_msg_list (list): A list of error logs collected from the initialization and execution of the workflow. This is a class attribute.
        workunits (List[WorkFlow]): A list of WorkUnit instance carried by this workflow.
        intermediate_data_mapper (dict): Runtime intermediate data for the current layer of workflow.
        inherited_data_mapper (dict): Data inherited from parent workflow.
        indexes (List[int]): A list of the hierarchical indexes of the work units to which the current workflow belongs in the configuration file.
                                  If the current workflow is at the outermost level, this value is an empty list.
        execution_status (int): A flag indicating the execution status of the WorkFlow instance. Using value from `class ExecutionStatus`.
    """
    debug: bool                 # Indicates whether to run in debug mode.
    error_msg_list = list()     # Error information. This is a class attribute.

    # The list of Work Unit.
    workunits: list
    
    # Intermediate Data.
    intermediate_data_mapper: dict  # Runtime intermediate data for the current layer of workflow.
    inherited_data_mapper: dict     # Data inherited from parent workflow.

    indexes: list  # A list of the hierarchical indexes of the work units to which the current workflow belongs in the configuration file.
    status: int   # A flag indicating the execution status of the WorkFlow instance.

    # Output
    # data_output_path = './Mutation.dat'

    @property
    def identifier(self):
        """The unique name of the WorkFlow instance, collapsed from `api_key` and `indexes`."""
        name = f"workflow@{'.'.join(map(str, self.indexes))}"
        return name
    
    @property
    def is_executable(self) -> bool:
        """Flag indicating whether the current instance passed self-inspection."""
        if self.status in StatusCode.unexecutable_statuses:
            return False
        else:
            return True
    
    #region WorkFlow Initialization.
    def __init__(self, debug: bool = False, 
                indexes: List[int] = list(), 
                intermediate_data_mapper: Dict[str, Any] = dict(), 
                inherited_data_mapper: Dict[str, Any] = dict()
                ) -> None:
        """Initialize an instance.
        
        Args:
            debug (bool, optional): Indicates whether to run in debug mode.
            unit_indexes (List[int]): A list of the hierarchical indexes of the work units to which the current workflow belongs in the configuration file.
                                      If the current workflow is at the outermost level, this value is an empty list.
            intermediate_data_mapper (Dict[str, Any], optional): The initial value of intermediate_data_mapper.
            inherited_data_mapper (Dict[str, Any], optional): Data Mapper containing data to be inherited.
        """
        # Assigned values.
        self.debug = debug
        self.indexes = indexes
        self.intermediate_data_mapper = intermediate_data_mapper
        self.inherited_data_mapper = inherited_data_mapper
        _LOGGER.debug(f"inherited_data_mapper: {self.inherited_data_mapper}")

        # Default values.
        self.workunits = list()
        self.status = StatusCode.CREATED
        return
    
    @classmethod
    def from_list(cls, units_dict_list: List[dict], debug: bool = False, 
                intermediate_data_mapper: Dict[str, Any] = dict(), 
                inherited_data_mapper: Dict[str, Any] = dict(), 
                indexes: List[int] = list()
                ) -> WorkFlow:
        """Initialize an instance of WorkFlow from a given dictionary instance.
        
        Args:
            units_dict_list (List[dict]): A list containing the WorkUnit configuration dicts.
            debug (bool, optional): Indicates whether to run in debug mode.
            unit_indexes (List[int]): A list of the hierarchical indexes of the work units to which the current workflow belongs in the configuration file.
                                      If the current workflow is at the outermost level, this value is an empty list.
            intermediate_data_mapper (Dict[str, Any], optional): The initial value of intermediate_data_mapper.
            inherited_data_mapper (Dict[str, Any], optional): Data Mapper containing data to be inherited.
        
        Returns:
            An instance of WorkFlow.
        """
        flow = cls(debug=debug, indexes=indexes, intermediate_data_mapper=intermediate_data_mapper, inherited_data_mapper=inherited_data_mapper)
        flow.status = StatusCode.INITIALIZING
        if (units_dict_list == None):
            no_workunit_err_msg = "Initializing WorkFlow with no workunits. Workunits are expected."
            flow.error_msg_list.append(no_workunit_err_msg)
            flow.status = StatusCode.FAILED_INITIALIZATION
        else:
            for unit_index, unit_dict in enumerate(units_dict_list):
                unit_indexes_to_pass = indexes.copy()
                unit_indexes_to_pass.append(unit_index)
                flow.add_unit(unit_dict, unit_indexes_to_pass)
                continue
        
        # If the flow status keeps INITIALIZING, then mark it as READY_TO_START.
        if flow.status == StatusCode.INITIALIZING:
            flow.status = StatusCode.READY_TO_START
        return flow
    
    @classmethod
    def from_json_string(cls, json_str: str, debug: bool = False) -> WorkFlow:
        """Initialize an instance of WorkFlow from a serialized json string.
        
        Args:
            json_str (str): The serialized json string.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        units_dict_list = loads(json_str)
        return WorkFlow.from_list(units_dict_list, debug=debug)

    @classmethod
    def from_json_file_object(cls, json_fobj: TextIOWrapper, debug: bool = False) -> WorkFlow:
        """Initialize an instance of WorkFlow from a json file object.
        
        Args:
            json_fobj (TextIOWrapper): The instance of a json file containing workflow configuration information.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        units_dict_list = load(json_fobj)
        return WorkFlow.from_list(units_dict_list, debug=debug)

    @classmethod
    def from_json_filepath(cls, json_filepath: str, debug: bool = False) -> WorkFlow:
        """Initialize an instance of WorkFlow from a json filepath.
        
        Args:
            json_filepath (str): The path to a json file containing workflow configuration information.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        with open(json_filepath) as fobj:
            return WorkFlow.from_json_file_object(fobj, debug=debug)
    
    def add_unit(self, unit_dict: dict, unit_indexes_to_pass: List[int] = list()) -> None:
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
            unit_indexes_to_pass (List[int]): A list of the hierarchical indexes of the workunit to add.
        
        Returns:
            The WorkUnit instance added.
        """
        workunit = WorkUnit(unit_dict, self, unit_indexes_to_pass, self.debug)   # Placeholder.
        api_key = workunit.api_key
        if (workunit.api_key in CONTROL_API_KEYS):
            _LOGGER.info(f"Parsing Control WorkUnit {workunit.identifier}.")
            if api_key == LOOP_API_KEY:
                # Setup a placeholder in `workflow.intermediate_data_mapper`.
                workunit = LoopWorkUnit.from_dict(unit_dict=unit_dict, workflow=self, indexes=unit_indexes_to_pass, debug=self.debug)
            elif api_key == GENERAL_API_KEY:
                error_msg = "GeneralWorkUnit should not be configured inside a workflow."
                _LOGGER.error(error_msg)
                self.error_msg_list.append(error_msg)
        else:
            _LOGGER.info(f"Parsing Science WorkUnit {workunit.identifier}.")
            workunit = WorkUnit.from_dict(unit_dict=unit_dict, workflow=self, indexes=unit_indexes_to_pass, debug=self.debug)

        # Setup a placeholder in `workflow.intermediate_data_mapper`.
        self.intermediate_data_mapper[VarFormatter.unformat_variable(workunit.return_key)] = PLACEHOLDER_VALUE

        # Add new workunit to the workflow.
        self.workunits.append(workunit)
        if workunit.status == StatusCode.FAILED_INITIALIZATION:
            self.status = StatusCode.FAILED_INITIALIZATION
            self.error_msg_list += workunit.error_msg_list
        return
    
    #endregion

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
            _LOGGER.error(FAILED_INITIALIZATION_ERROR_MSG)
            raise ValueError(self)
        else:
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

class WorkUnit():
    """Represents a single unit in a procedural workflow.

    - This class encapsulates the information and functionality for initializing, interpreting and executing
    a single unit of work in a procedural workflow. 
    - It handles API mappings, argument processing, and execution
    of a specified function or method.

    Attributes:
        indexes (list): A list of hierarchical indexes of the workunit.
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
        error_info_list (list): List of error messages encountered during processing.
        params_to_assign_at_execution (list): List of parameters that need to be assigned values at execution time.
        execution_status (int): A flag indicating the execution status of the WorkUnit instance.
        debug (bool): Flag indicating whether debugging mode is enabled.
    """
    indexes: list
    api_key: str
    args_dict_input: dict

    workflow: WorkFlow
    api: Callable
    return_key: str
    return_value: Any
    args_dict_to_pass: dict
    error_msg_list: list
    params_to_assign_at_execution: list
    status: int

    debug: bool

    @property
    def identifier(self) -> str:
        """The unique name of the WorkUnit instance, collapsed from `api_key` and `indexes`."""
        name = f"{self.api_key}@{':'.join(map(str, self.indexes))}"
        return name
    
    @property
    def is_executable(self) -> bool:
        """Flag indicating whether the current instance passed self-inspection."""
        if self.status in StatusCode.unexecutable_statuses:
            return False
        else:
            return True
        

    def __init__(self, unit_dict: dict = dict(), workflow: WorkFlow = None, indexes: List[int] = list(), debug: bool = False):
        """Initializes a WorkUnit instance.

        Args:
            unit_dict (dict): A dictionary derived from a JSON configuration, containing keys
                              for API mapping and execution.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            indexes (list, optional): A list of hierarchical indexes of the workunit.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        # Initialize attributes in the constructor
        # to prevent them from being considered as class attributes.
        self.api_key = unit_dict.get(WORKUNIT_API_NAME_KEY, str())
        self.return_key = unit_dict.get(WORKUNIT_RETURN_VALUE_KEY, None)
        self.args_dict_input = unit_dict.get(WORKUNIT_ARGUMENT_LIST_KEY, dict())
        self.args_dict_to_pass = dict()
        self.error_msg_list = list()
        self.params_to_assign_at_execution = list()
        self.return_value = None
        self.status = StatusCode.CREATED

        if workflow == None:
            self.workflow = WorkFlow()
        else:
            self.workflow = workflow
        self.indexes = indexes
        self.debug = debug
        return
    
    @classmethod
    def from_dict(cls, unit_dict: dict, workflow: WorkFlow = None, indexes: List[int] = list(), debug: bool = False) -> WorkUnit:
        """Initializes an instance of WorkUnit from a given dictionary.

        The dictionary should contain keys such as `api`, `store_as`, and `args`
        specifying the API to be called, the return key for storing the result, and the arguments
        for the API call, respectively.

        Args:
            unit_dict (dict): A dictionary derived from a JSON configuration, containing keys
                              for API mapping and execution.
            workflow (WorkFlow, optional): The workflow instance to which this unit belongs.
            indexes (list, optional): A list of hierarchical indexes of the workunit.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            WorkUnit: An initialized instance of WorkUnit.
        
        Raises:
            KeyError: If the API key cannot be mapped to any API.
        """
        # Initialize the class.
        unit = cls(unit_dict=unit_dict, workflow=workflow, indexes=indexes, debug=debug)
        unit.status = StatusCode.INITIALIZING

        # Map the API.
        unit.api = SCIENCE_API_MAPPER.get(unit.api_key, None)
        if unit.api is None:    # Inspect API Mapping.
            error_msg = f"Input API Key `{unit.api_key}` cannot be mapped to any API."
            _LOGGER.error(error_msg)
            if unit.debug:
                raise KeyError(unit.api_key)
            else:
                unit.error_msg_list.append(error_msg)
                return unit
        
        # Self Inspection.
        unit.self_inspection_and_args_reassembly()
        unit.check_initialization_status()

        return unit

    def execute(self) -> tuple:
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
            if self.return_key:
                self.return_value = self.api(**self.args_dict_to_pass)
            else:   # If no `return_key`, then `return_value` is None.
                self.api(**self.args_dict_to_pass)
                self.return_value = None
            self.status = StatusCode.EXIT_OK
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

    def update_args_value_from_data_mapper(self) -> None:
        """Update the `args_dict_to_pass` by mapping data from `workflow.intermediate_data_mapper` or `workflow.inherited_data_mapper`.
        If a value is mapped in both `inherited_data_mapper` and `intermediate_data_mapper`,
        then the latter one will overwrite the former one.
        """
        for param in self.params_to_assign_at_execution:
            for mapper in [self.workflow.intermediate_data_mapper, self.workflow.inherited_data_mapper]:
                if (mapped_datum := mapper.get(VarFormatter.unformat_variable(self.args_dict_input[param]))):
                    self.args_dict_to_pass[param] = mapped_datum
                    break
            continue
    
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

        # Checks if param.annotation is a typing type.
        is_typing_type = hasattr(annotation, '__origin__')

        if is_typing_type:
            # Check the origin of typing types (e.g., list for List[int])
            return isinstance(arg_value, annotation.__origin__)
        else:
            if (not isinstance(annotation, type)):
                # TODO Sometimes, the annotation may be parsed as a string value,
                # so we use pydoc.locate to cast it from string to type if it happens.
                index_of_bracket = annotation.find("[")
                matched_type = annotation[:index_of_bracket].lower() if index_of_bracket != -1 else annotation.lower()
                from pydoc import locate
                matched_type = locate(matched_type)

                if matched_type:
                    annotation = matched_type
                else:
                    annotation = eval(annotation)

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
                if (mapped_value):   # Check if there is a non-none-or-empty value in Mapper.
                    if (isinstance(mapped_value, type(PLACEHOLDER_VALUE)) # Placeholder values need to be excluded.
                        and mapped_value == PLACEHOLDER_VALUE):
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
        in `error_info_list`, and the `is_self_inspection_passed` flag is set to False.

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
            # if (WorkUnit.__inspect_argument_type(value, dict)):
            #     # `Dict` or `dict` type data are not inheritable.
            #     continue
            data_mapper_to_inherite[key] = value
        return data_mapper_to_inherite

class ControlWorkUnit(WorkUnit):
    """This class is an abstract class defining shared attributes and methods of Control WorkUnits (e.g. loop, general, etc.)
    This class is derived from WorkUnit class.
    """

    def __init__(self, unit_dict: Dict = dict(), workflow: WorkFlow = None, indexes: List[int] = list(), debug: bool = False):
        super().__init__(unit_dict, workflow, indexes, debug)

    def encode_layer_index(self, execution_index: int = None) -> str:
        """The index of current layer to be added to `indexes`.
        
        Args:
            execution_index (int, optional): Indicates the execution index of the target sub-workflow instance in the current control unit, 
                                            such as the index of loops in the `LoopWorkUnit`,
                                            and/or the index of parallel tasks in the `ParallelWorkUnit`. (TODO Zhong)
        
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
            layer_index (str): The index of current layer in `indexes`.

        Return:
            The execution index obtained. Default 0.
        """
        splitter = "_"
        if (splitter in layer_index):
            index_str = layer_index.split(splitter)[-1]
            index = int(index_str)
            return index
        else:
            return 0


class LoopWorkUnit(ControlWorkUnit):
    """This class is used to represent a special kind of workunit in a procedural workflow that is used to carry a loop.
    In this kind of workunit, all single units contained in the loop body is packaged together as a subworkflow.
    This class is derived from WorkUnit class.

    - This class encapsulates a loop structure within a workflow, where each iteration of the loop is treated as a sub-workflow.
    - It initializes and executes these sub-workflows successively. Failures in individual sub-workflows do not halt the entire loop's execution.

    Attributes:
        iterable_data (Iterable): The data over which the loop iterates.
        sub_workflow_list (List[WorkFlow]): A list of sub-workflows, each representing one iteration of the loop.
        Other inherited attributes...

    Methods:
        loop_unit_placeholder_api: A static placeholder API for initialization.
        __init__: Initializes a LoopWorkUnit instance.
        from_dict: Creates a LoopWorkUnit instance from a dictionary configuration.
        execute: Executes the loop, iterating over each sub-workflow.
    """

    iterable_data: Iterable
    sub_workflow_list: List[WorkFlow]

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

    def __init__(self, unit_dict: dict, workflow: WorkFlow = None, indexes: List[int] = list(), debug: bool = False):
        """Initializes a LoopWorkUnit instance.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            indexes (list, optional): A list of hierarchical indexes of the workunit.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        super().__init__(unit_dict, workflow, indexes, debug)
        self.sub_workflow_list = list()

    @classmethod
    def from_dict(cls, unit_dict: dict, workflow: WorkFlow = None, indexes: List[int] = list(), debug: bool = False) -> LoopWorkUnit:
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
            indexes (list, optional): A list of hierarchical indexes of the workunit.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            LoopWorkUnit: An instance of LoopWorkUnit initialized with the given configuration.
        """
        # Initialize the class.
        unit = cls(unit_dict, workflow, indexes, debug)
        unit.api = LoopWorkUnit.loop_unit_placeholder_api
        unit.status = StatusCode.INITIALIZING

        # Self Inspection.
        unit.self_inspection_and_args_reassembly()
        unit.check_initialization_status()

        # Initialize the Loop Body.
        _LOGGER.debug("Initializing a LoopWorkUnit now...")
        data_mapper_to_inherite = unit.generate_data_mapper_to_inherite()

        # Some data should be passed before the initialization of a sub-WorkFlow.
        loop_body_datum_varname = VarFormatter.unformat_variable(unit.args_dict_to_pass.get(LOOP_BODY_DATUM_LABEL))
        data_mapper_for_initialize_as_intermediate = {loop_body_datum_varname: PLACEHOLDER_VALUE}

        # Initialize a PlaceHolder Sub-WorkFlow to perform self-inspection.
        flow_indexes = unit.indexes.copy()
        flow_indexes.append(unit.encode_layer_index())
        placeholder_sub_workflow = WorkFlow.from_list(units_dict_list=unit.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL),
            debug=unit.debug, intermediate_data_mapper=data_mapper_for_initialize_as_intermediate,
            inherited_data_mapper=data_mapper_to_inherite, indexes=flow_indexes)
        
        # Add error message from placeholder workflow initialization.
        if (placeholder_sub_workflow.error_msg_list):
            unit.error_msg_list += placeholder_sub_workflow.error_msg_list

        return unit

    def execute(self) -> tuple:
        """Executes the loop, iterating over each sub-workflow in `sub_workflow_list`.

        In each iteration, the corresponding sub-workflow is executed, and its results are collected.
        Failures in sub-workflows do not stop the execution of the entire loop.

        Returns:
            tuple: A tuple containing the return key and a list of results from each loop iteration.
        """
        self.return_value = list()
        data_mapper_to_inherite = self.generate_data_mapper_to_inherite()
        self.update_args_value_from_data_mapper()
        _LOGGER.info(f'We are about to execute a LoopWorkUnit ...')
        self.status = StatusCode.RUNNING
        loop_body_datum_varname = self.args_dict_to_pass.get(LOOP_BODY_DATUM_LABEL)
        for loop_index, loop_datum in enumerate(self.args_dict_to_pass.get(LOOP_ITERABLE_DATA_LABEL)):
            
            # Initialize the sub-workflow at execution time.
            data_mapper_for_initialization = {loop_body_datum_varname: loop_datum}
            flow_indexes = self.indexes.copy()
            flow_indexes.append(self.encode_layer_index(loop_index))
            workflow = WorkFlow.from_list(
                units_dict_list=self.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL),
                debug=self.debug, intermediate_data_mapper=data_mapper_for_initialization,
                inherited_data_mapper=data_mapper_to_inherite, indexes=flow_indexes)

            # Execute the sub-workflow.
            _LOGGER.info(f'Executing loop_{loop_index} ...')
            self.sub_workflow_list.append(workflow)
            workflow_return = workflow.execute()
            self.return_value.append({self.encode_layer_index(loop_index): workflow_return})

            self.update_execution_status(workflow=workflow)
            continue

        if (self.status == StatusCode.RUNNING):
            self.status = StatusCode.EXIT_OK
        return self.return_key, self.return_value
    
    def update_execution_status(self, workflow: WorkFlow) -> None:
        """Update the execution status of LoopWorkUnit instance.
        
        Args:
            workflow (WorkFlow): The WorkFlow instance from loop bodies.
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
                self.status = StatusCode.EXIT_WITH_ERROR
            return
        elif workflow.status in StatusCode.pause_excluding_error_statuses:
            if self.status in StatusCode.error_including_pause_statuses:
                self.status = StatusCode.EXIT_WITH_ERROR_AND_PAUSE
            else:
                self.status = StatusCode.RUNNING_WITH_PAUSE_IN_INNER_UNITS
            return

class GeneralWorkUnit(ControlWorkUnit):
    """This class is used to represent a special kind of unit of work that would normally be used at the outermost level 
    of a workflow configuration, existing as a carrier for the entire workflow.
    This class is derived from WorkUnit class.

    Attributes:
        sub_workflow (WorkFlow): The workflow it carries.
        Other inherited attributes...

    Methods:
        general_unit_placeholder_api: A static placeholder API for initialization.
        __init__: Initializes a GeneralWorkUnit instance.
        from_dict: Creates a GeneralWorkUnit instance from a dictionary configuration.
        execute: Executes the GeneralWorkUnit, i.e., executes the sub_workflow.
    """
    sub_workflow: WorkFlow
    working_directory: str
    autodump: bool

    @staticmethod
    def general_unit_placeholder_api(workunits: list, **kwargs):
        """Serves as a placeholder API for initializing LoopWorkUnit. This API is not intended for actual execution but marks
        the necessary arguments for initialization.

        Args:
            workunits (list): A list containing the WorkUnit configuration dicts.
        """
        return

    def __init__(self, unit_dict: dict, workflow: WorkFlow = None, working_directory: str = DEFAULT_WORK_DIR, autodump: bool = False, debug: bool = False):
        """Initializes a GeneralWorkUnit instance.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            debug (bool, optional): Indicates whether to run in debug mode.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to. Usually None.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            autodump (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
        """
        super().__init__(unit_dict=unit_dict, workflow=workflow, debug=debug)
        self.working_directory = working_directory
        self.autodump = autodump

    @classmethod
    def from_dict(cls, unit_dict: dict, working_directory: str = DEFAULT_WORK_DIR, autodump: bool = False, debug: bool = False) -> GeneralWorkUnit:
        """Initializes an instance of GeneralWorkUnit from a given dictionary.

        Here, however, `api` must be `general` to indicate that the WorkUnit is a GeneralWorkUnit, 
        and the arguments contained in `args` are an assemble of WorkUnits forming a WorkFlow and 
        required information to execute the GeneralWorkUnit (the whole workflow).

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            autodump (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            LoopWorkUnit: An instance of LoopWorkUnit initialized with the given configuration.
        """
        # Initialize the class.
        unit = cls(unit_dict=unit_dict, working_directory=working_directory, autodump=autodump, debug=debug)
        if not autodump:
            not_autodump_warning_msg = "The `autodump` value is not set to `True`, so this GeneralWorkUnit instance will not automatically save its state on error or pause and continue computing later after corrections."
            _LOGGER.warning(not_autodump_warning_msg)
            _LOGGER.warning("Continue in 5 seconds...\n")
            sleep(5.0)

        unit.api = GeneralWorkUnit.general_unit_placeholder_api

        # Self Inspection.
        unit.self_inspection_and_args_reassembly()
        unit.check_initialization_status()
        
        # Initialize the WorkFlow.
        data_mapper_for_init = unit.args_dict_to_pass.copy()
        del data_mapper_for_init[CONTROL_BODY_WORKUNITS_LABEL]  # Pass kwargs to the workflow. Kwargs contain information specifying CPU/GPU partitions.

        if (units_dict_list:=unit.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL)):
            unit.sub_workflow = WorkFlow.from_list(units_dict_list, debug=debug, intermediate_data_mapper=data_mapper_for_init)
        else:
            error_msg = "Initializing GeneralWorkUnit with empty `workunits`. WorkUnit(s) are expected."
            _LOGGER.error(error_msg)
            raise ValueError(error_msg)
        return unit
    
    @classmethod
    def from_json_string(cls, json_str: str, working_directory: str = DEFAULT_WORK_DIR, autodump: bool = False, debug: bool = False) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a serialized json string.
        
        Args:
            json_str (str): The serialized json string.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            autodump (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            debug (bool, optional): Indicates whether to run in debug mode. Default False.
        
        Returns:
            An instance of WorkFlow.
        """
        unit_dict = loads(json_str)
        return GeneralWorkUnit.from_dict(unit_dict=unit_dict, working_directory=working_directory, autodump=autodump, debug=debug)

    @classmethod
    def from_json_file_object(cls, json_fobj: TextIOWrapper, working_directory: str = DEFAULT_WORK_DIR, autodump: bool = False, debug: bool = False) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a json file object.
        
        Args:
            json_fobj (TextIOWrapper): The instance of a json file containing GeneralWorkUnit configuration information.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            autodump (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            debug (bool, optional): Indicates whether to run in debug mode. Default False.
        
        Returns:
            An instance of GeneralWorkUnit.
        """
        unit_dict = load(json_fobj)
        return GeneralWorkUnit.from_dict(unit_dict=unit_dict, working_directory=working_directory, autodump=autodump, debug=debug)

    @classmethod
    def from_json_filepath(cls, json_filepath: str, working_directory: str = DEFAULT_WORK_DIR, autodump: bool = False, debug: bool = False) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a json filepath.
        
        Args:
            json_filepath (str): The path to a json file containing GeneralWorkUnit configuration information.
            working_directory (str, optional): Indicates where to run the job. Default to current working directory.
            autodump (bool, optional): Whether to automatically save the GeneralWorkUnit instance as a pickle file 
                                        when it exits due to Error, Pause, or Completion. Default False.
            debug (bool, optional): Indicates whether to run in debug mode. Default False.
        
        Returns:
            An instance of WorkFlow.
        """
        with open(json_filepath) as fobj:
            return GeneralWorkUnit.from_json_file_object(fobj, working_directory=working_directory, autodump=autodump, debug=debug)
        
    def execute(self) -> tuple:
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
            self.return_value = self.sub_workflow.execute()
            if (self.sub_workflow.status in StatusCode.error_excluding_pause_statuses):
                self.status = StatusCode.EXIT_WITH_ERROR_IN_INNER_UNITS
            elif ((self.sub_workflow.status in StatusCode.pause_excluding_error_statuses) and (self.status not in StatusCode.error_excluding_pause_statuses)):
                self.status = StatusCode.RUNNING_WITH_PAUSE_IN_INNER_UNITS
            else:
                self.status = StatusCode.EXIT_OK
            if (self.autodump):
                pickle_filepath = self.dump_snapshot_file()
                _LOGGER.info(f"The current running status has been saved to '{pickle_filepath}'. Please go to this location to view your file(s).")
            return self.return_key, self.return_value
    
    def dump_snapshot_file(self, save_directory: str = str()) -> str:
        """
        Saves the current state of the GeneralWorkUnit instance to a snapshot file.

        This method serializes the entire state of the GeneralWorkUnit instance and
        saves it to a file in the specified directory. The filename is generated
        based on the current timestamp.

        Note:
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
                pickle.dump(self, file=fobj, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as exc:
            if dir_to_save != self.working_directory:
                _LOGGER.warning(f"Unable to save file to '{dir_to_save}'. Using your working directory '{self.working_directory}' instead.")
                dir_to_save = self.working_directory
                filepath = path.join(dir_to_save, filename)
                with open(file=filepath, mode="wb") as fobj:
                    pickle.dump(self, file=fobj, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                raise exc
        return filepath
    
    @staticmethod
    def load_snapshot_file(filepath: str) -> GeneralWorkUnit:
        """
        Loads a GeneralWorkUnit instance from a snapshot file.

        This method deserializes the state of a GeneralWorkUnit instance from the given
        snapshot file. If the deserialized object is not an instance of GeneralWorkUnit,
        it raises a TypeError.

        Args:
            filepath (str): The path to the snapshot file from which to load the state.

        Returns:
            GeneralWorkUnit: The deserialized GeneralWorkUnit instance.

        Raises:
            TypeError: If the deserialized object is not an instance of GeneralWorkUnit.
        """
        with open(file=filepath, mode="rb") as fobj:
            unit = pickle.load(file=fobj)
            
            if not isinstance(unit, GeneralWorkUnit):
                raise TypeError("The loaded object is not an instance of GeneralWorkUnit.")
            
            return unit