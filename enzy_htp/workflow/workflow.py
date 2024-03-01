#! python3
# -*- encoding: utf-8 -*-
"""
WorkFlow and WorkUnit.

TODO (Zhong): Typing type annotation inspection. (Completed 2024-02-16 23:16 UTC-6)
TODO (Zhong): Annotation inspection in intermediate_data_mapper. (Completed 2024-03-01 01:00 UTC-6)
TODO (Zhong): APIs entitled `loop`. (Completed 2024-02-26 20:01 UTC-6)
TODO (Zhong): How to handle intermediate data transfer between parent and child workflows. Child workflow should be able to read the data in parent flow. (Completed 2024-02-28 06:06 UTC-6)
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
from json import load, loads
from typing import Any, Dict, List, Type, Union
from collections.abc import Iterable, Callable

from enzy_htp.core.logger import _LOGGER

from .config import SCIENCE_API_MAPPER, JsonConfigVariableFormatter as VarFormatter

# API Keys for flow control, such as `loop`.
CONTROL_API_KEYS = [
    "loop",
    "general",
]

# In the Json file used to define the workflow, 
# keys in the mapper are defined by the developer,
# which are not user-definable.
WORKUNIT_API_NAME_KEY = "api"               # Key indicating the API name.
WORKUNIT_RETURN_VALUE_KEY = "store_as"      # Key indicating the key where the return value of the workunit is stored.
WORKUNIT_ARGUMENT_LIST_KEY = "args"         # Key indicates the arguments to be passed into the API.
CONTROL_BODY_WORKUNITS_LABEL = "workunits"  # Key indicating the list of workunits contained in the body in a unit playing a control role.
LOOP_ITERABLE_DATA_LABEL = "loop_data"      # Key indicating the object to iterate over in a LoopWorkUnit.
LOOP_BODY_DATUM_LABEL = "loop_datum"        # Key indicating the element iterated from ITERABLE_DATA in each workunit of the loop body.

PLACEHOLDER_VALUE = "PLAC$H0LD$R"   # Value for placeholder.

class WorkFlow():
    """Rrepresents a procedural WorkFlow primarily consists of a number of WorkUnit instances, 
    including data mappers, self-inspection flags, to ensure its operation.

    Each WorkFlow Instance will perform self-inspection during initialization. If the self-inspection fails, it cannot be executed.
    
    If an error occurs during the execution of a single WorkUnit instance, the workflow will record the error, 
    jump out of the WorkUnit instance, and execute the next one until the execution of the entire workflow is completed.
    
    Attributes:
        debug (bool): Flag indicating whether debugging mode is enabled.
        error_msg_list (list): A list of error logs collected from the initialization and execution of the workflow.
        is_self_inspection_passed (bool): Indicates if the workflow has passed self inspection. If False, the workflow cannot be executed.
        workunits (List[WorkFlow]): A list of WorkUnit instance carried by this workflow.
        intermediate_data_mapper (dict): Runtime intermediate data for the current layer of workflow.
        inherited_data_mapper (dict): Data inherited from parent workflow.
    """
    debug: bool                 # Indicates whether to run in debug mode.
    error_msg_list = list()     # Error information.
    is_self_inspection_passed: bool

    # The list of Work Unit.
    workunits: list
    
    # Intermediate Data.
    intermediate_data_mapper: dict  # Runtime intermediate data for the current layer of workflow.
    inherited_data_mapper: dict     # Data inherited from parent workflow.

    # Output
    # data_output_path = './Mutation.dat'

    #region WorkFlow Initialization.
    def __init__(self, debug: bool = False, intermediate_data_mapper: Dict[str, Any] = dict(), inherited_data_mapper: Dict[str, Any] = dict()) -> None:
        """Initialize an instance.
        
        Args:
            debug (bool, optional): Indicates whether to run in debug mode.
            intermediate_data_mapper (Dict[str, Any], optional): The initial value of intermediate_data_mapper.
            inherited_data_mapper (Dict[str, Any], optional): Data Mapper containing data to be inherited.
        """
        self.debug = debug
        self.intermediate_data_mapper = intermediate_data_mapper
        self.inherited_data_mapper = inherited_data_mapper
        _LOGGER.debug(f"inherited_data_mapper: {self.inherited_data_mapper}")

        self.workunits = list()
        self.is_self_inspection_passed = False
        return
    
    @classmethod
    def from_list(cls, units_dict_list: List[dict], debug: bool = False, intermediate_data_mapper: Dict[str, Any] = dict(), inherited_data_mapper: Dict[str, Any] = dict()) -> WorkFlow:
        """Initialize an instance of WorkFlow from a given dictionary instance.
        
        Args:
            units_dict_list (List[dict]): A list containing the WorkUnit configuration dicts.
            debug (bool, optional): Indicates whether to run in debug mode.
            intermediate_data_mapper (Dict[str, Any], optional): The initial value of intermediate_data_mapper.
            inherited_data_mapper (Dict[str, Any], optional): Data Mapper containing data to be inherited.
        
        Returns:
            An instance of WorkFlow.
        """
        flow = cls(debug=debug, intermediate_data_mapper=intermediate_data_mapper, inherited_data_mapper=inherited_data_mapper)
        if (units_dict_list == None):
            flow.error_msg_list.append('Initializing WorkFlow with no workunits. Workunits are expected.')
        else:
            for unit_dict in units_dict_list:
                flow.add_unit(unit_dict)
                continue
        
        # If there's no error message, then self-inspection is passed.
        if (len(flow.error_msg_list) == 0):
            flow.is_self_inspection_passed = True
        else:
            pass
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
        api_key = unit_dict.get(WORKUNIT_API_NAME_KEY, str())
        workunit = WorkUnit()   # Placeholder.
        if (api_key in CONTROL_API_KEYS):
            _LOGGER.info(f"Parsing new Control WorkUnit: {api_key}.")
            if api_key == "loop":
                # Setup a placeholder in `workflow.intermediate_data_mapper`.
                workunit = LoopWorkUnit.from_dict(unit_dict=unit_dict, placeholder_workflow=self, debug=self.debug)
            elif api_key == "general":
                error_msg = "GeneralWorkUnit should not be configured inside a workflow."
                _LOGGER.error(error_msg)
                self.error_msg_list.append(error_msg)
        else:
            _LOGGER.info(f"Parsing new Science WorkUnit: {api_key}.")
            workunit = WorkUnit.from_dict(unit_dict=unit_dict, workflow=self, debug=self.debug)

        # Setup a placeholder in `workflow.intermediate_data_mapper`.
        self.intermediate_data_mapper[VarFormatter.unformat_variable(workunit.return_key)] = PLACEHOLDER_VALUE

        # Add new workunit to the workflow.
        self.workunits.append(workunit)
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
        if not self.is_self_inspection_passed:
            _LOGGER.error('The initialized WorkFlow has not yet passed the self-inspection, so it is not allowed to be executed.')
            raise ValueError(self)
        else:
            _LOGGER.info('The initialized WorkFlow has successfully passed the self-inspection. Proceeding to execution.')
            for workunit in self.workunits:
                workunit: WorkUnit
                _LOGGER.info(f'Executing WorkUnit API Key: {workunit.api_key} ...')
                return_key, return_value = workunit.execute()
                self.intermediate_data_mapper[VarFormatter.unformat_variable(return_key)] = return_value
                continue
            return self.intermediate_data_mapper


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
    api_key: str
    args_dict_input: dict

    workflow: WorkFlow
    api: Callable
    return_key: str
    return_value: Any
    args_dict_to_pass: dict
    error_msg_list: list
    is_self_inspection_passed: bool
    params_to_assign_at_execution: list

    debug: bool

    def __init__(self, workflow: WorkFlow = None, debug: bool = False):
        """Initializes a WorkUnit instance.

        Args:
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        # Initialize attributes in the constructor
        # to prevent them from being considered as class attributes.
        self.api_key = str()
        self.args_dict_input = dict()
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
        unit.api_key = unit_dict.get(WORKUNIT_API_NAME_KEY, str())
        unit.return_key = unit_dict.get(WORKUNIT_RETURN_VALUE_KEY, None)
        unit.args_dict_input = unit_dict.get(WORKUNIT_ARGUMENT_LIST_KEY, dict())

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
        unit.check_self_inspection_result()

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
            return self.return_key, self.return_value
        except Exception as exc:
            _LOGGER.error(f'Catching Error when executing `{self.api_key}`:')
            if (self.debug):
                raise exc
            else:
                _LOGGER.error(exc)
                self.error_msg_list.append(str(exc))
                return self.return_key, None

    def update_args_value_from_data_mapper(self) -> None:
        """Update the `args_dict_to_pass` by mapping data from `intermediate_data_mapper` or `inherited_data_mapper`.
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
                for key, value in self.args_dict_input.items():
                    if key not in self.args_dict_to_pass.keys():
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
                    is_inspection_passed = False
                else:   # Skip if it's optional.
                    pass
                continue

            # Map data from `intermediate_data_mapper` or `inherited_data_mapper`.
            _LOGGER.debug(f"{self.api_key}, param_name: {param.name}, arg_value: {arg_value}")
            if VarFormatter.is_variable_formatted(arg_value):   # Check if the arg_value is a formatted as a value name.
                unformatted_arg_value = VarFormatter.unformat_variable(arg_value)
                is_key_mapped = False   # Indicate if the value is mapped.
                is_value_mapped_none = True # Indicate if the value in the mapper is None at present.
                for mapper in [self.workflow.intermediate_data_mapper, self.workflow.inherited_data_mapper]:    # Listed in order of priority.
                    if (mapper and unformatted_arg_value in mapper):    # If mapped, leave for later.
                        is_key_mapped = True
                        self.params_to_assign_at_execution.append(param_name)
                        self.args_dict_to_pass[param.name] = unformatted_arg_value
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
                    continue
                elif (not is_key_mapped):
                    error_msg = f"The param `{param_name}` in API `{self.api_key} hasn't received expected variable {arg_value}`. Check your configuration and there's maybe some typos!"
                    _LOGGER.error(error_msg)
                    self.error_msg_list.append(error_msg)
                    continue
                _LOGGER.debug(f"{self.api_key}, param_name: {param.name}, arg_value: {unformatted_arg_value}, Value Changed.")
            
            # Inspect argument type.
            if (param.annotation != _empty):
                annotation = param.annotation
                pass_type_inspection = WorkUnit.__inspect_argument_type(arg_value, annotation)
                if not pass_type_inspection:
                    error_msg = f'Receiving argument `{param.name}` with value `{arg_value}` in unexpected type {type(arg_value)} while {annotation} is expected (when initializing `{self.api_key}`).'
                    _LOGGER.error(error_msg)
                    unit_error_list.append(error_msg)
                    is_inspection_passed = False

            # Inspection passed, now loading arguments.
            self.args_dict_to_pass[param.name] = arg_value
            continue
        
        self.is_self_inspection_passed = is_inspection_passed
        self.error_msg_list += unit_error_list
        return
    
    def check_self_inspection_result(self) -> None:
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
                raise ValueError(self)
            else:
                pass
        else:
            return

    def generate_data_mapper_to_inherite(self) -> Dict[str, Any]:
        """Generate data mappers for sub-workflow inheritance.
        
        `Dict` or `dict` type data are not inheritable.
        
        Returns:
            Data mapper to inherite (dict)
        """
        data_mapper_to_inherite = dict()
        if (inherited_data_mapper:=self.workflow.inherited_data_mapper):
            # Data Mapper inherited from parent workflow(s) should be inherited.
            data_mapper_to_inherite.update(inherited_data_mapper.copy())
        for key, value in self.workflow.intermediate_data_mapper.items():
            if (WorkUnit.__inspect_argument_type(value, dict)):
                # `Dict` or `dict` type data are not inheritable.
                continue
            data_mapper_to_inherite[key] = value
        return data_mapper_to_inherite

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
    def loop_unit_placeholder_api(loop_data: Iterable, workunits: list):
        """Serves as a placeholder API for initializing LoopWorkUnit. This API is not intended for actual execution but marks
        the necessary arguments for initialization.

        Args:
            loop_data (Iterable): The data to be iterated over in the loop.
            workunits (dict): A list containing the WorkUnit configuration dicts.
        """
        pass

    def __init__(self, workflow: WorkFlow = None, debug: bool = False):
        """Initializes a LoopWorkUnit instance.

        Args:
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            debug (bool, optional): Indicates whether to run in debug mode.
        """
        super().__init__(workflow, debug)
        self.sub_workflow_list = list()

    @classmethod
    def from_dict(cls, unit_dict: dict, placeholder_workflow: WorkFlow = None, debug: bool = False) -> LoopWorkUnit:
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
            debug (bool, optional): Indicates whether to run in debug mode.

        Returns:
            LoopWorkUnit: An instance of LoopWorkUnit initialized with the given configuration.
        """
        # Initialize the class.
        unit = cls(placeholder_workflow, debug=debug)
        unit.api_key = unit_dict.get(WORKUNIT_API_NAME_KEY, str())
        unit.return_key = unit_dict.get(WORKUNIT_RETURN_VALUE_KEY, None)
        unit.args_dict_input = unit_dict.get(WORKUNIT_ARGUMENT_LIST_KEY, dict())
        unit.api = LoopWorkUnit.loop_unit_placeholder_api

        # Self Inspection.
        unit.self_inspection_and_args_reassembly()
        unit.check_self_inspection_result()

        # Initialize the Loop Body.
        _LOGGER.debug("Initializing a LoopWorkUnit now...")
        data_mapper_to_inherite = unit.generate_data_mapper_to_inherite()

        # Some data should be passed before the initialization of a sub-WorkFlow.
        data_mapper_for_initialization = {LOOP_BODY_DATUM_LABEL: PLACEHOLDER_VALUE}

        # Initialize a PlaceHolder Sub-WorkFlow to perform self-inspection.
        placeholder_workflow = WorkFlow.from_list(units_dict_list=unit.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL),
            debug=unit.debug, intermediate_data_mapper=data_mapper_for_initialization,
            inherited_data_mapper=data_mapper_to_inherite)
        
        # Add error message from placeholder workflow initialization.
        if (placeholder_workflow.error_msg_list):
            unit.error_msg_list += placeholder_workflow.error_msg_list

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
        data_mapper_to_inherite = self.generate_data_mapper_to_inherite()
        self.update_args_value_from_data_mapper()
        _LOGGER.info(f'We are about to execute a LoopWorkUnit ...')
        for loop_datum in self.args_dict_to_pass.get(LOOP_ITERABLE_DATA_LABEL):
            
            # Initialize the sub-workflow at execution time.
            data_mapper_for_initialization = {LOOP_BODY_DATUM_LABEL: loop_datum}
            workflow = WorkFlow.from_list(
                units_dict_list=self.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL),
                debug=self.debug, intermediate_data_mapper=data_mapper_for_initialization,
                inherited_data_mapper=data_mapper_to_inherite)

            # Execute the sub-workflow.
            _LOGGER.info(f'Executing loop_{loop_index} ...')
            workflow_return = workflow.execute()
            self.sub_workflow_list.append(workflow)
            self.return_value.append({f'loop_{loop_index}': workflow_return})
            loop_index += 1
        return self.return_key, self.return_value

class GeneralWorkUnit(WorkUnit):
    """This class is used to represent a special kind of unit of work that would normally be used at the outermost level 
    of a workflow configuration, existing as a carrier for the entire workflow.
    This class is derived from WorkUnit class.

    Attributes:
        workflow (List[WorkFlow]): The workflow it carries.
        Other inherited attributes...

    Methods:
        general_unit_placeholder_api: A static placeholder API for initialization.
        __init__: Initializes a GeneralWorkUnit instance.
        from_dict: Creates a GeneralWorkUnit instance from a dictionary configuration.
        execute: Executes the GeneralWorkUnit, i.e., executes the workflow.
    """
    workflow: WorkFlow

    @staticmethod
    def general_unit_placeholder_api(workunits: list, **kwargs):
        """Serves as a placeholder API for initializing LoopWorkUnit. This API is not intended for actual execution but marks
        the necessary arguments for initialization.

        Args:
            workunits (list): A list containing the WorkUnit configuration dicts.
        """
        return

    def __init__(self, debug: bool = False, workflow: WorkFlow = None):
        """Initializes a GeneralWorkUnit instance.

        Args:
            debug (bool, optional): Indicates whether to run in debug mode.
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
        """
        super().__init__(workflow, debug)

    @classmethod
    def from_dict(cls, unit_dict: dict, debug: bool = False, workflow: WorkFlow = None) -> GeneralWorkUnit:
        """Initializes an instance of GeneralWorkUnit from a given dictionary.

        Here, however, `api` must be `general` to indicate that the WorkUnit is a GeneralWorkUnit, 
        and the arguments contained in `args` are an assemble of WorkUnits forming a WorkFlow and 
        required information to execute the GeneralWorkUnit (the whole workflow).

        Args:
            unit_dict (dict): Configuration dictionary with necessary parameters for initialization.
            debug (bool, optional): Indicates whether to run in debug mode.
            workflow (WorkFlow, optional): The parent workflow of this unit. Usually None.

        Returns:
            LoopWorkUnit: An instance of LoopWorkUnit initialized with the given configuration.
        """
        # Initialize the class.
        unit = cls(debug=debug, workflow=workflow)
        unit.api_key = unit_dict.get(WORKUNIT_API_NAME_KEY, str())
        unit.return_key = unit_dict.get(WORKUNIT_RETURN_VALUE_KEY, None)
        unit.args_dict_input = unit_dict.get(WORKUNIT_ARGUMENT_LIST_KEY, dict())
        unit.api = GeneralWorkUnit.general_unit_placeholder_api

        # Self Inspection.
        unit.self_inspection_and_args_reassembly()
        unit.check_self_inspection_result()
        
        # Initialize the WorkFlow.
        data_mapper_for_init = unit.args_dict_to_pass.copy()
        del data_mapper_for_init[CONTROL_BODY_WORKUNITS_LABEL]  # Pass kwargs to the workflow. Kwargs contain information specifying CPU/GPU partitions.

        if (units_dict_list:=unit.args_dict_to_pass.get(CONTROL_BODY_WORKUNITS_LABEL)):
            unit.workflow = WorkFlow.from_list(units_dict_list, debug=debug, intermediate_data_mapper=data_mapper_for_init)
        else:
            error_msg = "Initializing GeneralWorkUnit with empty `workunits`. WorkUnit(s) are expected."
            _LOGGER.error(error_msg)
            raise ValueError(error_msg)
        return unit
    
    @classmethod
    def from_json_string(cls, json_str: str, debug: bool = False) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a serialized json string.
        
        Args:
            json_str (str): The serialized json string.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        unit_dict = loads(json_str)
        return WorkFlow.from_list(unit_dict, debug=debug)

    @classmethod
    def from_json_file_object(cls, json_fobj: TextIOWrapper, debug: bool = False) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a json file object.
        
        Args:
            json_fobj (TextIOWrapper): The instance of a json file containing GeneralWorkUnit configuration information.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of GeneralWorkUnit.
        """
        unit_dict = load(json_fobj)
        return GeneralWorkUnit.from_dict(unit_dict, debug=debug)

    @classmethod
    def from_json_filepath(cls, json_filepath: str, debug: bool = False) -> GeneralWorkUnit:
        """Initialize an instance of GeneralWorkUnit from a json filepath.
        
        Args:
            json_filepath (str): The path to a json file containing GeneralWorkUnit configuration information.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        """
        with open(json_filepath) as fobj:
            return GeneralWorkUnit.from_json_file_object(fobj, debug=debug)
        
    def execute(self) -> tuple:
        """Executes the GeneralWorkUnit.

        Returns:
            tuple: A tuple containing the return key and the intermediate data from the workflow.
        """
        self.return_value = self.workflow.execute()
        return self.return_key, self.return_value
