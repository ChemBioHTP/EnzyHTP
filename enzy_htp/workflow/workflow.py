#! python3
# -*- encoding: utf-8 -*-
"""
WorkFlow and WorkUnit.

TODO (Zhong): Typing type annotation inspection. (Completed 2024-02-16 23:16 UTC-6)
TODO (Zhong): Annotation inspection in intermediate_data_mapper.
TODO (Zhong): APIs entitled `loop`. (Completed 2024-02-26 20:01 UTC-6)
TODO (Zhong): How to handle intermediate data transfer between parent and child workflows.
TODO (Zhong): Data Output.

@File    :   workflow.py
@Created :   2024/02/05 22:53
@Author  :   Zhong, Yinjie
@Contact :   yinjie.zhong@vanderbilt.edu
"""

# Here put the import lib.
from __future__ import annotations
from io import TextIOWrapper, FileIO
from inspect import _empty, Parameter, signature
from json import load
from typing import Any, Dict, List, Type, Union
from collections.abc import Iterable, Callable

from enzy_htp.core.logger import _LOGGER

from . import SCIENCE_API_MAPPER

# API Keys for flow control, such as `loop`.
CONTROL_API_KEYS = [
    "loop"
]

# In the Json file used to define the workflow, 
# keys in the mapper are defined by the developer,
# which are not user-definable.
WORKUNIT_API_NAME_KEY = "api"           # Key indicating the API name.
WORKUNIT_RETURN_VALUE_KEY = "store_as"  # Key indicating the key where the return value of the workunit is stored.
WORKUNIT_ARGUMENT_LIST_KEY = "args"     # Key indicates the arguments to be passed into the API.
LOOP_ITERABLE_DATA_LABEL = "iterable_data"  # Key indicating the object to iterate over in a LoopWorkUnit.
LOOP_BODY_WORKUNITS_LABEL = "actions"       # Key indicating the list of workunits contained in the loop body in a LoopWorkUnit.
LOOP_BODY_DATUM_LABEL = "loop_body_datum"   # Key indicating the element iterated from ITERABLE_DATA in each workunit of the loop body.

class WorkFlow():
    """Class to interpret/compose the workflow from/to a json file.
    
    Attributes:
        wildtype_pdb_path (str): The path to a wildtype pdb file (Default: Empty).
        cpu_account (str): The account with access to CPU nodes.
        cpu_partition (str): The CPU Partition to use.
        gpu_account (str): The account with access to GPU nodes.
        gpu_partition (str): The GPU Partition to use.
    """
    execution_type: str = str()
    parallel_groups: int = 1
    mut_parallel_groups: int = 1
    debug: bool = False

    # Error information.
    error_msg_list = list()
    is_self_inspection_passed: bool = False

    # The list of Work Unit.
    workunits = list()
    
    # Intermediate Data.
    intermediate_data_mapper: dict

    # Output
    # data_output_path = './Mutation.dat'

    #region WorkFlow Initialization.
    def __init__(self, debug: bool = False, intermediate_data_mapper: Dict[str, Any] = dict()) -> None:
        """Initialize an instance.
        
        Args:
            debug (bool, optional): Indicates whether to run in debug mode.
            intermediate_data_mapper (Dict[str, Any], optional): The initial value of intermediate_data_mapper.
        """
        self.debug = debug
        self.intermediate_data_mapper = intermediate_data_mapper
        return
    
    @classmethod
    def from_dict(cls, flow_dict: dict, debug: bool = False, intermediate_data_mapper: Dict[str, Any] = dict()) -> WorkFlow:
        """Initialize an instance of WorkUnit from a given dictionary instance.
        
        Args:
            flow_dict (dict): A dictionary probably derived from config json file.
            debug (bool, optional): Indicates whether to run in debug mode.
            intermediate_data_mapper (Dict[str, Any], optional): The initial value of intermediate_data_mapper.
        
        Returns:
            An instance of WorkFlow.
        """
        flow = cls(debug=debug, intermediate_data_mapper=intermediate_data_mapper)
        flow.execution_type = flow_dict.get('execution_type', str())
        flow.parallel_groups = flow_dict.get('parallel_groups', 1)
        flow.mut_parallel_groups = flow_dict.get('mut_parallel_groups', 1)
        
        unit_dict_list = flow_dict.get('workunits', None)
        if (unit_dict_list == None):
            flow.error_msg_list.append('Initializing WorkFlow with no workunits. Workunits are expected.')
        else:
            for unit_dict in unit_dict_list:
                flow.add_unit(unit_dict)
                continue
        
        # If there's no error message, then self-inspection is passed.
        if (len(flow.error_msg_list) == 0):
            flow.is_self_inspection_passed = True
        else:
            [_LOGGER.error(error_msg) for error_msg in flow.error_msg_list]
        return flow

    @classmethod
    def from_json_file_object(cls, json_fobj: TextIOWrapper, debug: bool = False) -> WorkFlow:
        """Initialize an instance of WorkUnit from a given dictionary instance.
        
        Args:
            json_fobj (TextIOWrapper): The instance of a json file containing workflow configuration information.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        flow_dict = load(json_fobj)
        return WorkFlow.from_dict(flow_dict, debug=debug)

    @classmethod
    def from_json_filepath(cls, json_filepath: str, debug: bool = False) -> WorkFlow:
        """Initialize an instance of WorkUnit from a json filepath.
        
        Args:
            json_filepath (str): The path to a json file containing workflow configuration information.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        with open(json_filepath) as fobj:
            return WorkFlow.from_json_file_object(fobj, debug=debug)
    
    def add_unit(self, unit_dict: dict) -> None:
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
        """
        if (api_key:=unit_dict.get(WORKUNIT_API_NAME_KEY, str()) in CONTROL_API_KEYS):
            if api_key == "loop":
                workunit = LoopWorkUnit.from_dict(unit_dict=unit_dict, workflow=self, debug=self.debug)
        else:
            workunit = WorkUnit.from_dict(unit_dict=unit_dict, workflow=self, debug=self.debug)
            self.error_msg_list += workunit.error_msg_list
            self.workunits.append(workunit)
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
        if not self.is_self_inspection_passed:
            _LOGGER.error('The initialized WorkFlow has not yet passed the self-inspection, so it is not allowed to be executed.')
            raise ValueError(self)
        else:
            _LOGGER.info('The initialized WorkFlow has successfully passed the self-inspection. Proceeding to execution.')
            for workunit in self.workunits:
                workunit: WorkUnit
                return_key, return_value = workunit.execute()
                self.intermediate_data_mapper[return_key] = return_value
                continue


class WorkUnit():
    """Represents a single unit in a procedural workflow.

    - This class encapsulates the information and functionality for initializing, interpreting and executing
    a single unit of work in a procedural workflow. 
    - It handles API mappings, argument processing, and execution
    of a specified function or method.

    Attributes:
        __api_key (str): Key used to identify and map the appropriate science API function.
                         Read from the `api` field in the configuration.
        __args_dict_input (dict): The raw, unprocessed arguments dictionary as specified
                                  in each unit of the JSON configuration.
        workflow (WorkFlow): Reference to the parent workflow object that this unit is part of.
        api (object): The actual science API function object mapped using `__api_key`.
        return_key (str): The key under which the return value of the API call will be stored.
                          Read from the `store_as` field in the configuration.
        return_value: The value returned by the execution of the API function.
        args_dict_to_pass (dict): Processed arguments dictionary to be passed to the API function.
        error_info_list (list): List of error messages encountered during processing.
        is_self_inspection_passed (bool): Flag indicating whether self-inspection passed.
        params_to_assign_at_execution (list): List of parameters that need to be assigned values at execution time.
        debug (bool): Flag indicating whether debugging mode is enabled.
    """
    __api_key: str
    __args_dict_input: dict

    workflow: WorkFlow
    api: Callable
    return_key: str
    return_value: Any
    args_dict_to_pass: dict
    error_msg_list: list
    is_self_inspection_passed: bool
    params_to_assign_at_execution: list

    debug: bool = False

    def __init__(self, workflow: WorkFlow = None, debug: bool = False):
        """Initializes a WorkUnit instance.

        Args:
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        # Initialize attributes in the constructor
        # to prevent them from being considered as class attributes.
        self.__api_key = str()
        self.__args_dict_input = dict()
        self.args_dict_to_pass = dict()
        self.error_msg_list = list()
        self.is_self_inspection_passed = False
        self.params_to_assign_at_execution = list()
        self.return_value = None

        if workflow == None:
            self.workflow = WorkFlow()
        else:
            self.workflow = workflow
        self.debug = debug
        return
    
    @classmethod
    def from_dict(cls, unit_dict: dict, workflow: WorkFlow = None, debug: bool = False) -> WorkUnit:
        """Initializes an instance of WorkUnit from a given dictionary.

        The dictionary should contain keys such as `api`, `store_as`, and `args`
        specifying the API to be called, the return key for storing the result, and the arguments
        for the API call, respectively.

        Args:
            unit_dict (dict): A dictionary derived from a JSON configuration, containing keys
                              for API mapping and execution.
            workflow (WorkFlow, optional): The workflow instance to which this unit belongs.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            WorkUnit: An initialized instance of WorkUnit.
        
        Raises:
            KeyError: If the API key cannot be mapped to any API.
        """
        # Initialize the class.
        unit = cls(workflow, debug=debug)
        unit.__api_key = unit_dict.get(WORKUNIT_API_NAME_KEY, str())
        unit.return_key = unit_dict.get(WORKUNIT_RETURN_VALUE_KEY, None)
        unit.__args_dict_input = unit_dict.get(WORKUNIT_ARGUMENT_LIST_KEY, dict())

        # Map the API.
        unit.api = SCIENCE_API_MAPPER.get(unit.__api_key, None)
        if unit.api is None:    # Inspect API Mapping.
            error_msg = f"Input API Key `{unit.__api_key}` cannot be mapped to any API."
            if unit.debug:
                _LOGGER.error(error_msg)
                raise KeyError(unit.__api_key)
            else:
                unit.error_msg_list.append(error_msg)
                return unit
        
        # Self Inspection.
        unit.__self_inspection_and_reassembly()
        unit.__check_self_inspection_result()

        # Setup a placeholder in `workflow.intermediate_data_mapper`.
        unit.workflow.intermediate_data_mapper[unit.return_key] = None
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
            for param in self.params_to_assign_at_execution:
                self.args_dict_to_pass[param] = self.workflow.intermediate_data_mapper.get(self.__args_dict_input[param])
            if self.return_key:
                self.return_value = self.api(**self.args_dict_to_pass)
            else:   # If no `return_key`, then `return_value` is None.
                self.api(**self.args_dict_to_pass)
                self.return_value = None
            return self.return_key, self.return_value
        except Exception as exc:
            _LOGGER.error(f'Catching Error when executing `{self.__api_key}`:')
            if (self.debug):
                raise exc
            else:
                _LOGGER.error(exc)
                self.error_msg_list.append(str(exc))
                return self.return_key, None
    
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
            # Inspect basic types.
            return isinstance(arg_value, annotation)

    def __self_inspection_and_reassembly(self) -> None:
        """Performs self-inspection and reassembles the arguments for the API call.

        This method is called every time a new instance is created. It performs a self-inspection
        to validate the `args` dictionary specified in the JSON configuration. It checks for the
        existence of required arguments and their types. If any errors are found, they are recorded
        in `error_info_list`, and the `is_self_inspection_passed` flag is set to False.

        The arguments that pass the inspection are reassembled into `args_dict_to_pass` for the API call.
        """
        is_inspection_passed = True
        unit_error_list = list()
        
        # Inspect Arguments.
        api_parameters = signature(self.api).parameters
        for param_name in api_parameters:
            param: Parameter = api_parameters[param_name]

            # If `*args` or `**kwargs` parameters exist (which is always at the very last),
            # add all the arguments to the `args_dict_to_pass`. Then, break the loop.
            # APIs with both `*args` and `**kwargs` are not and will not be supported.
            if (param.name in ['arg', 'args', 'kwarg', 'kwargs']):
                for key, value in self.__args_dict_input.items():
                    if key not in self.args_dict_to_pass.keys():
                        self.args_dict_to_pass[key] = value
                        continue
                break

            # A param with empty default value is a required parameter.
            is_required = False
            if param.default == _empty:
                is_required = True

            # Inspect argument existence.
            arg_value = self.__args_dict_input.get(param.name, None)
            if not arg_value:   # Missing the argument.
                if is_required: # Record error if it's required.
                    unit_error_list.append(f'Missing required argument `{param.name}` in API `{self.__api_key}`.')
                    is_inspection_passed = False
                else:   # Skip if it's optional.
                    pass
                continue

            # Map intermediate data from `intermediate_data_mapper`.
            if isinstance(arg_value, str):
                if (arg_value in self.workflow.intermediate_data_mapper):  # If mapped, leave for later.
                    self.params_to_assign_at_execution.append(param_name)
                    self.args_dict_to_pass[param.name] = arg_value
                    continue
            
            # Inspect argument type.
            if (param.annotation != _empty):
                pass_type_inspection = WorkUnit.__inspect_argument_type(arg_value, param.annotation)
                if not pass_type_inspection:
                    unit_error_list.append(f'Receiving argument `{param.name}` in unexpected type {type(arg_value)} while {param.annotation} is expected (when initializing `{self.__api_key}`).')
                    is_inspection_passed = False

            # Inspection passed, now loading arguments.
            self.args_dict_to_pass[param.name] = arg_value
            continue
        
        self.is_self_inspection_passed = is_inspection_passed
        self.error_msg_list += unit_error_list
        return
    
    def __check_self_inspection_result(self) -> None:
        """Checks the result of self-inspection of the instance.

        This method checks if the self-inspection of the instance passed (`is_self_inspection_ok`). 
        If the inspection fails and debugging is enabled, it logs each error message and raises a 
        ValueError with the list of error messages. If debugging is not enabled, it appends the 
        error messages to the workflow error list.

        Raises:
            ValueError: If self-inspection fails and debugging is enabled, a ValueError is raised 
                        with all collected error information.
        """
        if self.is_self_inspection_passed == False:
            if self.debug:
                for error_msg in self.error_msg_list:
                    _LOGGER.error(error_msg)
                raise ValueError(self)
            else:
                pass
        else:
            return


class LoopWorkUnit(WorkUnit):
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
    def loop_unit_placeholder_api(iterable_data: Iterable, actions: dict):
        """Serves as a placeholder API for initializing LoopWorkUnit. This API is not intended for actual execution but marks
        the necessary arguments for initialization.

        Args:
            iterable_data (Iterable): The data to be iterated over in the loop.
            actions (dict): A dictionary defining the actions (sub-workflow) for each iteration.
        """
        pass

    def __init__(self, workflow: WorkFlow = None, debug: bool = False):
        """Initializes a LoopWorkUnit instance.

        Args:
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        super().__init__(workflow, debug)

    @classmethod
    def from_dict(cls, unit_dict: dict, workflow: WorkFlow = None, debug: bool = False) -> LoopWorkUnit:
        """Initializes an instance of LoopWorkUnit from a given dictionary.

        The dictionary should contain keys such as `api`, `store_as`, and `args`
        specifying the API to be called, the return key for storing the result, and the arguments
        for the API call, respectively.

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            workflow (WorkFlow, optional): The parent workflow of this unit.
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            LoopWorkUnit: An instance of LoopWorkUnit initialized with the given configuration.
        """
        # Initialize the class.
        unit = cls(workflow, debug=debug)
        unit.__api_key = unit_dict.get(WORKUNIT_API_NAME_KEY, str())
        unit.return_key = unit_dict.get(WORKUNIT_RETURN_VALUE_KEY, None)
        unit.__args_dict_input = unit_dict.get(WORKUNIT_ARGUMENT_LIST_KEY, dict())
        unit.api = LoopWorkUnit.loop_unit_placeholder_api

        # Self Inspection.
        unit.__self_inspection_and_reassembly()
        unit.__check_self_inspection_result()

        if (iterable_data := unit.args_dict_to_pass.get(LOOP_ITERABLE_DATA_LABEL)):
            for loop_datum in iterable_data:
                # Some data should be passed before the initialization of a sub-WorkFlow.
                data_mapper_for_init = {LOOP_BODY_DATUM_LABEL: loop_datum}

                # Initialize the SubWorkflow.
                workflow = WorkFlow.from_dict(flow_dict=unit.args_dict_to_pass.get(LOOP_BODY_WORKUNITS_LABEL),
                    debug=unit.debug, intermediate_data_mapper=data_mapper_for_init)
                
                # Add the initialized sub-workflow to `sub_workflow_list`.
                unit.sub_workflow_list.append(workflow)
        else:
            unit.error_msg_list.append('Initializing LoopWorkUnit with empty iterable data list. Iterable data are expected.')

        # Setup a placeholder in `workflow.intermediate_data_mapper`.
        unit.workflow.intermediate_data_mapper[unit.return_key] = None
        return unit

    def execute(self) -> tuple:
        """Executes the loop, iterating over each sub-workflow in `sub_workflow_list`.

        In each iteration, the corresponding sub-workflow is executed, and its results are collected.
        Failures in sub-workflows do not stop the execution of the entire loop.

        Returns:
            tuple: A tuple containing the return key and a list of results from each loop iteration.
        """
        loop_index = 0
        self.return_value = list()
        for workflow in self.sub_workflow_list:
            _LOGGER.info(f'Executing loop_{loop_index} ...')
            workflow: WorkFlow
            workflow.execute()
            self.return_value.append({f'loop_{loop_index}': workflow.intermediate_data_mapper})
            loop_index += 1
        return self.return_key, self.return_value
