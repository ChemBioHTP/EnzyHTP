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

from .workflow_config import SCIENCE_API_MAPPER

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
        api_key (str): API Key in SCIENCE_API_MAPPER.
        api (object): Science API mapped by `api_key`.
        return_value (str): The file name/path returned by the API.
        args_dict_input (dict): Uninspected `args` dict set in the each procedure unit of json.
        args_dict (dict): A dictionary of inspected arguments for the API will be loaded.
    '''
    api_key: str
    api: object
    return_value: str
    args_dict_input: dict
    args_dict: dict

    def __init__(self):
        '''Declare a unit of the procedure.'''
        self.args_dict = dict()
        pass
    
    @classmethod
    def from_dict(cls, unit_dict: dict) -> ProcedureUnit:
        '''Initialize an instance of ProcedureUnit from a given dictionary containing keys entitled `api`, `return` and `args`.
        
        Args:
            unit_dict (dict): A dictionary derived from the loaded json file/string, which should contain keys entitled `api`, `return` and `args`.
        
        Returns:
            An instance of ProcedureUnit.
        '''
        unit = cls()
        unit.api_key = unit_dict.get('api', str())
        unit.api = SCIENCE_API_MAPPER.get(unit.api_key, None)
        if unit.api is None:    # Inspect API Mapping.
            error_info = f'Input API Key `{unit.api_key}` cannot be mapped to any api.'
            _LOGGER.error(error_info)
            raise KeyError(error_info)
        
        unit.return_value = unit_dict.get('return', str())
        unit.args_dict_input = unit_dict.get('args', dict())
        
        is_args_correct, error_info_list = unit.self_inspection_and_reassembly()
        if (is_args_correct == False):
            for error_info in error_info_list:
                _LOGGER.error(error_info)
            raise ValueError(error_info_list)
        else:
            return unit

    def execute(self):
        '''Execute the procedure unit.'''
        return self.api(**self.args_dict)

    def self_inspection_and_reassembly(self) -> tuple(bool, list):
        '''Perform a self-inspection to check if the `args` dict set in json is correct.
        Then, reassemble the arguments to args_dict.
        
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
                elif (arg_value == None):   # Skip if it's optional.
                    pass
                continue

            if (param.annotation is not _empty) and (type(arg_value) is not param.annotation):    # Unexpected argument type.
                error_info_list.append(f'Receiving argument `{param.name}` in unexpected type {type(arg_value)} (while {param.annotation} is expected).')
                is_pass = False
                continue

            # TODO (Zhong): *args and **kwargs parameters should be handled.

            # Inspection passed, now loading arguments.
            self.args_dict[param.name] = arg_value
            continue
        
        return (is_pass, error_info_list)


if __name__ == '__main__':
    cpu_account = 'yang_lab_csb'
    cpu_partition = 'production'
    gpu_account = 'csb_gpu_acc'
    gpu_partition= 'a6000x4'
    gpu_walltime = '5-00:00:00'