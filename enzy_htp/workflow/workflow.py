#! python3
# -*- encoding: utf-8 -*-
'''
WorkFlow and WorkUnit.

TODO (Zhong): Use `getfullargspec(func).annotations.get('return').__args__` to extract return annotations from a function and then perform type inspection.
TODO (Zhong): APIs entitled `loop`.
TODO (Zhong): Data Output.

@File    :   workflow.py
@Created :   2024/02/05 22:53
@Author  :   Zhong, Yinjie
@Contact :   yinjie.zhong@vanderbilt.edu
'''

# Here put the import lib.
from __future__ import annotations
from io import TextIOWrapper, FileIO
from inspect import _empty, Parameter, signature
from json import load
from typing import Union, overload

from enzy_htp.core.logger import _LOGGER

from .config import SCIENCE_API_MAPPER

class WorkFlow():
    '''Class to interpret/compose the workflow from/to a json file.
    
    Attributes:
        wildtype_pdb_path (str): The path to a wildtype pdb file (Default: Empty).
        cpu_account (str): The account with access to CPU nodes.
        cpu_partition (str): The CPU Partition to use.
        gpu_account (str): The account with access to GPU nodes.
        gpu_partition (str): The GPU Partition to use.
    '''
    execution_type: str = str()
    parallel_groups: int = 1
    mut_parallel_groups: int = 1

    # cpu_account: str = str()
    # cpu_partition: str = str()
    # gpu_account: str = str()
    # gpu_partition: str = str()
    # gpu_walltime: str = str()

    debug: bool = False

    # Error information.
    error_msg_list = list()
    is_self_inspection_passed: bool = False

    # The list of Work Unit.
    workunits = list()
    
    # Intermediate Data.
    intermediate_data_mapper = dict()

    # Output
    data_output_path = './Mutation.dat'

    #region WorkFlow Initialization.
    def __init__(self, debug: bool = False) -> None:
        '''Initialize an instance.
        
        Args:
            debug (bool, optional): Indicates whether to run in debug mode.
        '''
        self.debug = debug
        return
    
    @classmethod
    def from_dict(cls, flow_dict: dict, debug: bool = False) -> WorkFlow:
        '''Initialize an instance of WorkUnit from a given dictionary instance.
        
        Args:
            flow_dict (dict): A dictionary probably derived from config json file.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        '''
        flow = cls(debug=debug)
        flow.execution_type = flow_dict.get('execution_type', str())
        flow.parallel_groups = flow_dict.get('parallel_groups', 1)
        flow.mut_parallel_groups = flow_dict.get('mut_parallel_groups', 1)
        
        unit_dict_list = flow_dict.get('workunits', None)
        if (unit_dict_list == None):
            flow.error_msg_list.append('Initializing WorkFlow with no workunits. workunits are expected.')
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
        '''Initialize an instance of WorkUnit from a given dictionary instance.
        
        Args:
            json_fobj (TextIOWrapper): The instance of a json file containing workflow configuration information.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        '''
        flow_dict = load(json_fobj)
        return WorkFlow.from_dict(flow_dict, debug=debug)

    @classmethod
    def from_json_filepath(cls, json_filepath: str, debug: bool = False) -> WorkFlow:
        '''Initialize an instance of WorkUnit from a json filepath.
        
        Args:
            json_filepath (str): The path to a json file containing workflow configuration information.
            debug (bool, optional): Indicates whether to run in debug mode.
        
        Returns:
            An instance of WorkFlow.
        '''
        with open(json_filepath) as fobj:
            return WorkFlow.from_json_file_object(fobj, debug=debug)
    
    def add_unit(self, unit_dict: dict) -> None:
        '''Adds a new unit to the workflow.

        This method creates a new instance of WorkUnit from a provided dictionary and adds it to 
        the workflow. It initializes the WorkUnit with the given configuration `unit_dict`, 
        links it to the current workflow instance, and sets the debug mode according to the 
        provided argument. Any error messages encountered during the initialization of the WorkUnit 
        are appended to the workflow's error message list.

        Args:
            unit_dict (dict): A dictionary containing the configuration for the unit to be added.
                            This dictionary should include keys for API mapping and execution 
                            parameters.
        '''
        workunit = WorkUnit.from_dict(unit_dict=unit_dict, workflow=self, debug=self.debug)
        self.error_msg_list += workunit.error_msg_list
        self.workunits.append(workunit)
        return
    
    #endregion

    def execute(self):
        '''Executes all units in the workflow sequentially.

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

        Raises:
            ValueError: If the workflow has not passed the self-inspection.
        '''
        if not self.is_self_inspection_passed:
            _LOGGER.error('The initialized WorkFlow has not yet passed the self-inspection, so it is not allowed to be executed.')
            raise ValueError(self)
        else:
            _LOGGER.info('The initialized WorkFlow has successfully passed the self-inspection. Proceeding to execution.')
            for workunit in self.workunits:
                workunit: WorkUnit
                return_key, return_value = workunit.execute()
                self.intermediate_data_mapper[return_key] = return_value
                self.error_msg_list += workunit.error_msg_list
                continue


class WorkUnit():
    '''Represents a single unit in a procedural workflow.

    - This class encapsulates the information and functionality for initializing, interpreting and executing
    a single unit of work in a larger workflow. 
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
    '''
    __api_key: str = str()
    __args_dict_input: dict = dict()

    workflow: WorkFlow
    api: object
    return_key: str
    return_value = None
    args_dict_to_pass: dict
    error_msg_list: list
    is_self_inspection_passed: bool
    params_to_assign_at_execution: list

    debug: bool = False


    def __init__(self, workflow: WorkFlow = None, debug: bool = False):
        '''Initializes a WorkUnit instance.

        Args:
            workflow (WorkFlow): A reference to the workflow object this unit belongs to.
            debug (bool, optional): Indicates whether to run in debug mode.
        '''
        # Initialize attributes in the constructor
        # to prevent them from being considered as class attributes.
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
        '''Initializes an instance of WorkUnit from a given dictionary.

        The dictionary should contain keys such as `api`, `return` (or `store_as`), and `args`
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
        '''
        # Initialize the class.
        unit = cls(workflow, debug=debug)
        unit.__api_key = unit_dict.get('api', str())
        unit.return_key = unit_dict.get('store_as', str())
        unit.__args_dict_input = unit_dict.get('args', dict())

        # Map the API.
        unit.api = SCIENCE_API_MAPPER.get(unit.__api_key, None)
        if unit.api is None:    # Inspect API Mapping.
            error_msg = f'Input API Key `{unit.__api_key}` cannot be mapped to any API.'
            if unit.debug:
                _LOGGER.error(error_msg)
                raise KeyError(unit.__api_key)
            else:
                unit.error_msg_list.append(error_msg)
                return unit
        
        # Self Inspection.
        unit.self_inspection_and_reassembly()
        unit.check_self_inspection_result()

        # Setup a placeholder in `workflow.workunits_data_mapper`.
        # TODO (Zhong): Annotation inspection.
        unit.workflow.intermediate_data_mapper[unit.return_key] = None
        return unit

    def execute(self) -> tuple:
        '''Executes the workunit.

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
        '''
        try:
            for param in self.params_to_assign_at_execution:
                self.args_dict_to_pass[param] = self.workflow.intermediate_data_mapper.get(self.__args_dict_input[param])
            self.return_value = self.api(**self.args_dict_to_pass)
            self.workflow.intermediate_data_mapper[self.return_key] = self.return_value
            return self.return_key, self.return_value
        except Exception as exc:
            _LOGGER.error(f'Catching Error when executing `{self.__api_key}`:')
            if (self.debug):
                raise exc
            else:
                _LOGGER.error(exc)
                self.error_msg_list.append(str(exc))
                return self.return_key, None

    def self_inspection_and_reassembly(self) -> None:
        '''Performs self-inspection and reassembles the arguments for the API call.

        This method is called every time a new instance is created. It performs a self-inspection
        to validate the `args` dictionary specified in the JSON configuration. It checks for the
        existence of required arguments and their types. If any errors are found, they are recorded
        in `error_info_list`, and the `is_self_inspection_passed` flag is set to False.

        The arguments that pass the inspection are reassembled into `args_dict_to_pass` for the API call.
        '''
        is_inspection_passed = True
        unit_error_list = list()
        
        # Inspect Arguments.
        api_parameters = signature(self.api).parameters
        for param_name in api_parameters:
            param: Parameter = api_parameters[param_name]

            # If `*args` or `**kwargs` parameters exist (which is always at the very last),
            # add all the arguments to the `args_dict_to_pass`. Then, break the loop.
            # APIs with both `*args` and `**kwargs` are not and will not be supported.
            if (param.name in ['args', 'kwargs']):
                for key, value in self.__args_dict_input.items():
                    if key not in self.args_dict_to_pass.keys():
                        self.args_dict_to_pass[key] = value
                        continue
                break

            # A param with empty default value is a required parameter.
            is_required = False
            if param.default == _empty:
                is_required = True

            # Inspect argument existence and type.
            arg_value = self.__args_dict_input.get(param.name, None)
            if (arg_value == None):     # Missing the argument.
                if is_required: # Record error if it's required.
                    unit_error_list.append(f'Missing required argument `{param.name}` in API `{self.__api_key}`.')
                    is_inspection_passed = False
                else:   # Skip if it's optional.
                    pass
                continue

            # Map intermediate data from `workunits_data_mapper`.
            if isinstance(arg_value, str):
                if (arg_value in self.workflow.intermediate_data_mapper):  # If mapped, leave for later. TODO (Zhong): Annotation inspection.
                    self.params_to_assign_at_execution.append(param_name)
                    self.args_dict_to_pass[param.name] = arg_value
                    continue
            
            # Handle unexpected argument type.
            if (param.annotation != _empty) and (type(arg_value) != param.annotation):
                unit_error_list.append(f'Receiving argument `{param.name}` in unexpected type {type(arg_value)} while {param.annotation} is expected (when initializing `{self.__api_key}`).')
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
                for error_msg in self.error_msg_list:
                    _LOGGER.error(error_msg)
                raise ValueError(self)
            else:
                pass
        else:
            return

