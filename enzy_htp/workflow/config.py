"""this module interpret workflow configurations from config files
The configuration of the workflow is equivalent to the main script 
of the workflow (or shrapnel script) in term of information. This 
distinguish it from the _config module that this module 
also configs for the procedure of the workflow.

Author: QZ Shao <shaoqz@icloud.com>; Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Date: 2024-01-25"""

from collections.abc import Iterable
from typing import Any, List, Union

from enzy_htp import core, structure, preparation, mutation, geometry, quantum, analysis, interface, _LOGGER

SCIENCE_API_MAPPER = {
    "read_pdb" : structure.PDBParser().get_structure,
    "remove_solvent" : preparation.remove_solvent,
    "remove_hydrogens" : preparation.remove_hydrogens,
    "protonate_stru" : preparation.protonate_stru,
    "assign_mutant" : mutation.assign_mutant,
    "mutate_stru" : mutation.mutate_stru,
    "build_md_parameterizer": interface.amber.build_md_parameterizer,
    "equi_md_sampling" : geometry.equi_md_sampling,
    "single_point" : quantum.single_point,
    "optimize" : quantum.optimize,
    "ele_field_strength_at_along" : analysis.ele_field_strength_at_along,
    "ele_field_strength_at" : analysis.ele_field_strength_at,
    "ele_stab_energy_of_bond" : analysis.ele_stab_energy_of_bond,
    "ele_stab_energy_of_dipole" : analysis.ele_stab_energy_of_dipole,
    "bond_dipole": analysis.bond_dipole,
    "write_data_to_excel": core.write_data_to_excel,
}

class WidgetAPI:
    """
    Provides static methods for simple data processing tasks such as
    retrieving items and attributes from various Python objects.
    This utility class enhances data access in a generalized manner.
    """
    
    @staticmethod
    def get_item(obj: Union[list, tuple, dict, str], index: Any = 0) -> Any:
        """
        Retrieves an item from a given object (list, tuple, dictionary, or string)
        using the provided index or key.

        Args:
            obj (Union[list, tuple, dict, str]): The object from which to retrieve the item.
            index (Any): The index or key for accessing the item. Default is 0.

        Returns:
            Any: The item located at the specified index/key.
            If the obj is a dict and the index is an integer, a key-value pair will be returned.

        Raises:
            KeyError: If the key is not found in a dictionary.
            IndexError: If the index is out of range for lists or tuples.
            TypeError: If the object type does not support indexing.
        """
        try:
            if (isinstance(obj, dict) and isinstance(index, int)):
                dict_pairs_list = list(obj.items())
                return dict_pairs_list[index]
            return obj[index]
        except (KeyError, IndexError, TypeError) as e:
            raise e  # Optionally, handle or log the error here as needed.

    @staticmethod
    def get_attr(obj: Any, attr_chain: List[str]) -> Any:
        """
        Retrieves a nested attribute from an object based on a chain of attribute names.

        Args:
            obj (Any): The object from which attributes are to be fetched.
            attr_chain (List[str]): A list of attribute names representing the path to
                                    the desired attribute (e.g., ['attr1', 'attr2'] will fetch obj.attr1.attr2).

        Returns:
            Any: The value of the nested attribute, or None if any attribute in the chain does not exist.

        Example:
            If you pass an object `example` with a structure where `example.part1.subpart2.final` exists,
            calling `get_attr(example, ['part1', 'subpart2', 'final'])` will return the value of `final`.
        """
        attr = obj
        for attr_name in attr_chain:
            if hasattr(attr, attr_name):
                attr = getattr(attr, attr_name)
            else:
                raise AttributeError(f"The object of type '{type(attr)}' has no attribute '{attr_name}'. ")
        return attr


WIDGET_API_MAPPER = {
    "get_item": WidgetAPI.get_item,
    "get_attr": WidgetAPI.get_attr,
}

SCIENCE_API_MAPPER.update(WIDGET_API_MAPPER)

PARAM_METHOD_MAPPER = {
    "amber" : interface.amber.build_md_parameterizer,
}

class JsonConfigVariableFormatter():
    """Handles the formatting of variable names within JSON configuration files.

    In the JSON configuration file, if a parameter value is intended to be a variable name 
    instead of a literal value, it should be enclosed within backticks (`).
    
    For example: if you want to pass a variable named "pdb_path", you should annotate
    it as follows:
    
    ```json
    "args": {
          "path": "`pdb_path`"
    }
    ```    
    """

    @staticmethod
    def is_variable_formatted(value_to_check: Any) -> bool:
        """Check if a string is formatted as a variable name.

        Args:
            value_to_check (Any): Value to check for variable formatting.

        Returns:
            bool: True if the value_to_check is a string and is formatted as a variable name, False otherwise.
        """
        if (not isinstance(value_to_check, str)):
            return False
        return len(value_to_check) >= 3 and value_to_check.startswith("`") and value_to_check.endswith("`")
    
    @staticmethod
    def format_as_variable(name: str) -> str:
        """Format a string as a variable name by adding backticks.

        Args:
            name (str): The string to format as a variable name.

        Returns:
            str: The formatted variable name.
        """
        if not JsonConfigVariableFormatter.is_variable_formatted(name):
            return f"`{name}`"
        return name
        
    @staticmethod
    def unformat_variable(name: str) -> str:
        """Remove the backticks from a string formatted as a variable name.

        Args:
            name (str): The variable name to remove formatting from.

        Returns:
            str: The unformatted variable name.
        """
        if (JsonConfigVariableFormatter.is_variable_formatted(name)):
            name = name[1:] if name[0] == "`" else name
            name = name[:-1] if name[-1] == "`" else name
        return name

class Placeholder():
    """This class is to handle affairs relevant to Placeholder values."""
    # Placeholder values.
    PLACEHOLDER_STR_VALUE = "The world is yours and ours, but in the final analysis it's you young people's."  # Value for placeholder.
    PLACEHOLDER_INT_VALUE = -9876
    PLACEHOLDER_FLOAT_VALUE = 2024.03
    PLACEHOLDER_ITERABLE_VALUE = ["111", "565", "32162", "333", "555", "65132"]
    PLACEHOLDER_DICT_VALUE = {
        "zzz": "abbc",
        "yyy": "cdde",
        "xxx": "effg",
    }

    placeholder_values = [
        PLACEHOLDER_DICT_VALUE, 
        PLACEHOLDER_FLOAT_VALUE, 
        PLACEHOLDER_INT_VALUE, 
        PLACEHOLDER_ITERABLE_VALUE, 
        PLACEHOLDER_STR_VALUE
    ]
    
    @staticmethod
    def is_none_or_placeholder(var):
        """Check if a value is None or a Placeholder value.
        
        Args:
            var: The variable to check.

        Returns:
            The check result.
        """
        if var is None:
            return True
        else:
            placeholder_value = __class__.assign_placeholder_value(var)
            return var == placeholder_value

    @staticmethod
    def assign_placeholder_value(var):
        """
        Assigns var placeholder value to 'var' based on its type.

        Args:
            var: The variable whose placeholder value needs to be assigned.

        Returns:
            The placeholder value corresponding to the type of 'var'.
        """
        if isinstance(var, str):
            return __class__.PLACEHOLDER_STR_VALUE
        elif isinstance(var, int):
            return __class__.PLACEHOLDER_INT_VALUE
        elif isinstance(var, float):
            return __class__.PLACEHOLDER_FLOAT_VALUE
        elif isinstance(var, dict):
            return __class__.PLACEHOLDER_DICT_VALUE
        elif hasattr(var, '__iter__'):  # Checks if 'a' is iterable
            return __class__.PLACEHOLDER_ITERABLE_VALUE
        else:
            None

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
        READY_TO_START (int): The status when a WorkUnit or WorkFlow instance has passed the self-inspection but hasn't yet been started. Value: -6
        READY_WITH_UPDATES (int): The status when a previously executed WorkUnit or WorkFlow instance have its input arguments changed
                                due to the influence from the update in the input values of the task during continue computing. Value: -5
        SUSPECIOUS_UPDATES (int): The status when a previously `EXIT_OK` WorkUnit or WorkFlow instance has detected but unsure argument updates
                                    during reloading. Value: -4
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
        DEPRECATED (int): Indicates that the workunit or workflow is deprecated.
                        Any workflows or workunits inside it should be marked with this status and deleted. Value: 8
        FAILED_INITIALIZATION (int): Indicates that the initialization of the workunit or workflow failed. Value: 9
    """
    CREATED = -9
    INITIALIZING = -8
    READY_TO_START = -6
    READY_WITH_UPDATES = -5
    SUSPECIOUS_UPDATES = -4             # For Reload time only.
    RUNNING = -3
    RUNNING_WITH_PAUSE_IN_INNER_UNITS = -2   # For WorkFlow and ControlWorkUnit only.
    EXPECTED_PAUSE = -1                 # For Basic WorkUnit only. If a unit is paused, set its outer layers as `RUNNING_WITH_PAUSE_IN_INNER_UNITS`.
    EXIT_OK = 0
    EXIT_WITH_ERROR_IN_INNER_UNITS = 1  # For WorkFlow and ControlWorkUnit only.
    EXIT_WITH_ERROR = 2                 # For Science API only.
    EXIT_WITH_ERROR_AND_PAUSE = 3       # For WorkFlow and ControlWorkUnit only.
    DEPRECATED = 8
    FAILED_INITIALIZATION = 9

    pause_excluding_error_statuses = [EXPECTED_PAUSE, RUNNING_WITH_PAUSE_IN_INNER_UNITS]    # For logical judgment only.
    pause_including_error_statuses = [EXPECTED_PAUSE, RUNNING_WITH_PAUSE_IN_INNER_UNITS, EXIT_WITH_ERROR_AND_PAUSE] # For logical judgment only.
    error_excluding_pause_statuses = [EXIT_WITH_ERROR, EXIT_WITH_ERROR_IN_INNER_UNITS, FAILED_INITIALIZATION]       # For logical judgment only.
    error_including_pause_statuses = [EXIT_WITH_ERROR, EXIT_WITH_ERROR_IN_INNER_UNITS, FAILED_INITIALIZATION, EXIT_WITH_ERROR_AND_PAUSE]   # For logical judgment only.
    error_or_pause_statuses = [EXIT_WITH_ERROR, EXIT_WITH_ERROR_IN_INNER_UNITS, EXPECTED_PAUSE, RUNNING_WITH_PAUSE_IN_INNER_UNITS, EXIT_WITH_ERROR_AND_PAUSE]   # For logical judgment only.
    unexecutable_statuses = [CREATED, INITIALIZING, DEPRECATED, FAILED_INITIALIZATION]      # For logical judgment only.
    unexecuted_statuses = [CREATED, INITIALIZING, READY_TO_START, FAILED_INITIALIZATION]    # For logical judgment only. Note its distinction from `unexecutable_statuses`.
    skippable_statuses = [EXIT_OK]  # For logical judgement only.

demo_json_dict = {
    "workunits" : [
        {
            "api" : "read_pdb",
            "store_as" : "read_pdb_0",
            "args" : {
                "path" : "xxx",
            }
        },
        {
            "api" : "remove_hydrogens",
            "store_as" : "remove_hydrogens_0",
            "args" : {
                "stru" : "read_pdb_0",
                "polypeptide_only" : True,
            }
        },
        {
            "api" : "protonate_stru",
            "store_as" : "protonate_stru_0",
            "args" : {
                "stru" : "remove_hydrogens_0",
                "ph" : 7.0,
                "protonate_ligand" : False,
                "engine" : "pdb2pqr",
            }
        },
        {
            "api" : "assign_mutant",
            "store_as" : "assign_mutant_0",
            "args" : {
                "stru" : "protonate_stru_0",
                "pattern" : "xxx",
                "chain_sync_list" : ["A", "B"],
                "chain_index_mapper" : {"A" : 0, "B" : 100},
            }
        },
        {
            "api" : "loop", # make a for loop when read this
            "store_as" : "loop_0",
            "args" : {
                "data" : "assign_mutant_0",
                "actions" : [
                    {   
                        "api" : "mutate_stru",
                        "store_as" : "mutate_stru_0",
                        "args" : {
                            "stru" : "protonate_stru",
                            "mutant" : "loop_data",
                            "engine" : "pymol",
                            }
                    },
                    {
                        "api" : "equi_md_sampling",
                        "store_as" : "equi_md_sampling_0",
                        "args" : {
                            "stru" : "mutate_stru_0",
                            "param_method" : {
                                "engine" : "amber",
                                "force_fields" : "ff14SB",
                            },
                            "prod_time" : 100.0,
                            "prod_constrain" : [
                                {
                                    "type" : "distance",
                                    "atom_1" : "B.254.H2",
                                    "atom_2" : "A.101.OE2",
                                    "target_value" : 2.4,
                                },
                                {
                                    "type" : "angle",
                                    "atom_1" : "B.254.CAE",
                                    "atom_2" : "B.254.H2",
                                    "atom_3" : "A.101.OE2",
                                    "target_value" : 180.0,
                                },
                            ],
                            "cluster_job_config" : {
                                "cluster" : "slurm",
                                "res_keywords" : {
                                    "account" : "xxx",
                                    "partition" : "xxx"
                                }
                            }
                        }
                    },
                    {
                        "api" : "loop",
                        "store_as" : "loop_0_0",
                        "args" : {
                            "data" : "equi_md_sampling_0",
                            "actions" : [      
                                {   
                                    "api" : "loop",
                                    "store_as" : "loop_0_0_0",
                                    "args" : {
                                        "data" : "loop_data",
                                        "actions" : [ # target metrics goes here
                                            {   
                                                "api" : "ele_field_strength_at_along",
                                                "store_as" : "ef_0",
                                                "args" : {
                                                    "stru" : "loop_data",
                                                    "p1" : "A.1.CA",
                                                    "p2" : "A.1.CB",
                                                }
                                            },
                                        ]
                                    }
                                },

                            ]
                        }
                    }
                ],
            }
        },
    ],
    "execution_type" : "run", # or deploy
    "parallel_groups" : 5, # this will parallelize in the 1st loop.
    "mut_parallel_groups" : 5, # temp use for this version that works like shrapnel
}
"""TODO this could just be a python dict in pickle etc. In this case objects can be stored"""
