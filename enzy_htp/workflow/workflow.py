#! python3
# -*- encoding: utf-8 -*-
'''
@File    :   workflow.py
@Created :   2024/02/05 22:53
@Author  :   Zhong, Yinjie
@Contact :   yinjie.zhong@vanderbilt.edu
'''

# Here put the import lib.
from __future__ import annotations
from plum import dispatch
from inspect import _empty, Parameter, signature

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
    wildtype_pdb_path: str = str()
    cpu_account: str = str()
    cpu_partition: str = str()
    gpu_account: str = str()
    gpu_partition: str = str()
    gpu_walltime: str = str()
    
    # Intermediate Data.
    procedure_unit_return_values = dict()


    # Output
    data_output_path_dat = './Mutation.dat'

    def __init__(self, wildtype_pdb_path: str = str()) -> None:
        '''Initialize an instance.
        
        Args:
            wildtype_pdb_path (str): The path to a wildtype pdb file (Default: Empty).
        '''
        self.wildtype_pdb_path = wildtype_pdb_path
        return

class ProcedureUnit():
    '''Class for each unit in the procedure.
    
    Attributes:
        api_key (str): API Key in SCIENCE_API_MAPPER, read from `api` field.
        api (object): Science API mapped by `api_key`.
        return_key (str): The name that the return value from the api to be stored as, read from `store_as` field.
        return_value (str): The return value from the api.
        args_dict_input (dict): Uninspected `args` dict set in the each procedure unit of json.
        args_dict_to_pass (dict): A dictionary of inspected arguments for the API will be loaded.
    '''
    api_key: str = str()
    api: object = object()
    return_key: str = str()
    return_value = None
    args_dict_input: dict = dict()
    args_dict_to_pass: dict = dict()

    def __init__(self):
        '''Declare a unit of the procedure.'''
        self.args_dict_to_pass = dict()
        pass
    
    @classmethod
    def from_dict(cls, unit_dict: dict) -> ProcedureUnit:
        '''Initialize an instance of ProcedureUnit from a given dictionary containing keys entitled `api`, `return` and `args`.
        
        Args:
            unit_dict (dict): A dictionary derived from the loaded json file/string, which should contain keys entitled `api`, `return` and `args`.
        
        Returns:
            An instance of ProcedureUnit.
        '''
        # Initialize the class.
        unit = cls()
        unit.api_key = unit_dict.get('api', str())
        unit.return_key = unit_dict.get('store_as', str())
        unit.args_dict_input = unit_dict.get('args', dict())

        # Map the API.
        unit.api = SCIENCE_API_MAPPER.get(unit.api_key, None)
        if unit.api is None:    # Inspect API Mapping.
            error_info = f'Input API Key `{unit.api_key}` cannot be mapped to any API.'
            _LOGGER.error(error_info)
            raise KeyError(error_info)
        
        # Self Inspection.
        is_args_correct, error_info_list = unit.self_inspection_and_reassembly()
        if (is_args_correct == False):
            for error_info in error_info_list:
                _LOGGER.error(error_info)
            raise ValueError(error_info_list)
        else:
            return unit

    def execute(self):
        '''Execute the procedure unit.'''
        self.return_value = self.api(**self.args_dict_to_pass)
        return self.return_key, self.return_value

    def self_inspection_and_reassembly(self) -> tuple:
        '''Called everytime when a new instance is created.
        - Perform a self-inspection to check if the `args` dict set in json is correct.
        - Then, reassemble the arguments to `args_dict`.
        
        Returns:
            A tuple containing two elements:
            1. bool: Does the unit pass the self-inspection?
            2. list: Empty if the self-inspection is passed; Error information list if the self-inspection is failed.
        '''
        is_pass = True
        error_info_list = list()
        
        # Inspect Arguments.
        api_parameters = signature(self.api).parameters
        for param_name in api_parameters:
            param: Parameter = api_parameters[param_name]

            # If `*args` or `**kwargs` parameters exist (which is always at the very last),
            # add all the arguments to the `args_dict_to_pass`. Then, break the loop.
            # APIs with both `*args` and `**kwargs` are not and will not be supported.
            if (param.name in ['args', 'kwargs']):
                for key, value in self.args_dict_input.items():
                    if key not in self.args_dict_to_pass.keys():
                        self.args_dict_to_pass[key] = value
                        continue
                break

            # A param with empty default value is a required parameter.
            is_required = False
            if param.default is _empty:
                is_required = True

            # Inspect argument existence and type.
            arg_value = self.args_dict_input.get(param.name, None)
            if (arg_value == None):     # Missing the argument.
                if is_required: # Record error if it's required.
                    error_info_list.append(f'Missing required argument `{param.name}` in API `{self.api_key}`.')
                    is_pass = False
                else:   # Skip if it's optional.
                    pass
                continue

            if (param.annotation is not _empty) and (type(arg_value) is not param.annotation):    # Unexpected argument type.
                error_info_list.append(f'Receiving argument `{param.name}` in unexpected type {type(arg_value)} (while {param.annotation} is expected).')
                is_pass = False
                continue

            # Inspection passed, now loading arguments.
            self.args_dict_to_pass[param.name] = arg_value
            continue
        
        return (is_pass, error_info_list)


if __name__ == '__main__':
    cpu_account = 'yang_lab_csb'
    cpu_partition = 'production'
    gpu_account = 'csb_gpu_acc'
    gpu_partition= 'a6000x4'
    gpu_walltime = '5-00:00:00'