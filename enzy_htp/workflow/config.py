"""this module interpret workflow configurations from config files
The configuration of the workflow is equivalent to the main script 
of the workflow (or shrapnel script) in term of information. This 
distinguish it from the _config module that this module 
also configs for the procedure of the workflow.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-01-25"""

from enzy_htp import structure, preparation, mutation, geometry, quantum, analysis, interface
from typing import Any

SCIENCE_API_MAPPER = {
    "read_pdb" : structure.PDBParser().get_structure,
    "remove_solvent" : preparation.remove_solvent,
    "remove_hydrogens" : preparation.remove_hydrogens,
    "protonate_stru" : preparation.protonate_stru,
    "assign_mutant" : mutation.assign_mutant,
    "mutate_stru" : mutation.mutate_stru,
    "equi_md_sampling" : geometry.equi_md_sampling,
    "single_point" : quantum.single_point,
    "optimize" : quantum.optimize,
    "ele_field_strength_at_along" : analysis.ele_field_strength_at_along,
    "ele_field_strength_at" : analysis.ele_field_strength_at,
    "ele_stab_energy_of_bond" : analysis.ele_stab_energy_of_bond,
    "ele_stab_energy_of_dipole" : analysis.ele_stab_energy_of_dipole,
    "bond_dipole": analysis.bond_dipole,
}

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
