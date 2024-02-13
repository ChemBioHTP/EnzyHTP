"""
mutation module for enzy_htp.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-15
"""

from .mutation import (
    Mutation,
    generate_from_mutation_flag,
    check_repeat_mutation,
    remove_repeat_mutation,
    get_mutant_name_tag,
    get_involved_mutation,
)

from .api import (
    assign_mutant,
    mutate_stru,
)
