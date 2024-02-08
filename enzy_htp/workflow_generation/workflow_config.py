"""this module interpret workflow configurations from config files
The configuration of the workflow is equivalent to the main script 
of the workflow (or shrapnel script) in term of information. This 
distinguish it from the enzy_htp._config module that this module 
also configs for the procedure of the workflow.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-01-25"""

import enzy_htp

SCIENCE_API_MAPPER = {
    "read_pdb" : enzy_htp.structure.PDBParser().get_structure,
    "remove_solvent" : enzy_htp.preparation.remove_solvent,
    "remove_hydrogens" : enzy_htp.preparation.remove_hydrogens,
    "protonate_stru" : enzy_htp.preparation.protonate_stru,
    "assign_mutant" : enzy_htp.mutation.assign_mutant,
    "mutate_stru" : enzy_htp.mutation.mutate_stru,
    "equi_md_sampling" : enzy_htp.mutation.equi_md_sampling,
    "single_point" : enzy_htp.quantum.single_point,
    "optimize" : enzy_htp.quantum.optimize,
    "ele_field_strength_at_along" : enzy_htp.analysis.ele_field_strength_at_along,
    "ele_field_strength_at" : enzy_htp.analysis.ele_field_strength_at,
    "ele_stab_energy_of_bond" : enzy_htp.analysis.ele_stab_energy_of_bond,
    "ele_stab_energy_of_dipole" : enzy_htp.analysis.ele_stab_energy_of_dipole,
    "bond_dipole" : enzy_htp.analysis.bond_dipole,
}

PARAM_METHOD_MAPPER = {
    "amber" : enzy_htp.interface.amber.build_md_parameterizer,
}

demo_json = {
    "procedure" : [
        {
        "api" : "read_pdb",
        "name" : "read_pdb_0",
        "args" : {
            "path" : "xxx",
            }
        },
        {
        "api" : "remove_hydrogens",
        "name" : "remove_hydrogens_0",
        "args" : {
            "stru" : "read_pdb_0",
            "polypeptide_only" : True,
            }
        },
        {
        "api" : "protonate_stru",
        "name" : "protonate_stru_0",
        "args" : {
            "stru" : "remove_hydrogens_0",
            "ph" : 7.0,
            "protonate_ligand" : False,
            "engine" : "pdb2pqr",
            }
        },
        {
        "api" : "assign_mutant",
        "name" : "assign_mutant_0",
        "args" : {
            "stru" : "protonate_stru_0",
            "pattern" : "xxx",
            "chain_sync_list" : ["A", "B"],
            "chain_index_mapper" : {"A" : 0, "B" : 100},
            }
        },
        {
        "api" : "loop", # make a for loop when read this
        "name" : "loop_0",
        "args" : {
            "data" : "assign_mutant_0",
            "actions" : [
                {   
                "api" : "mutate_stru",
                "name" : "mutate_stru_0",
                "args" : {
                    "stru" : "protonate_stru",
                    "mutant" : "loop_0",
                    "engine" : "pymol",
                    }
                },
                {
                "api" : "equi_md_sampling",
                "name" : "equi_md_sampling_0",
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
                "name" : "loop_0_0",
                "args" : {
                    "data" : "equi_md_sampling_0",
                    "actions" : [      
                        {   
                        "api" : "loop",
                        "name" : "loop_0_0_0",
                        "args" : {
                            "data" : "loop_0_0",
                            "actions" : [ # target metrics goes here
                                {   
                                "api" : "ele_field_strength_at_along",
                                "name" : "ef_0",
                                "args" : {
                                    "stru" : "loop_0_0_0",
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
