"""
Mutation class for enzy_htp.


Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-02-13
"""

from .mutation import (
    Mutation,
    generate_from_mutation_flag,
    generate_mutation_from_target_list,
    check_repeat_mutation,
    remove_repeat_mutation,
    get_mutant_name_tag,
    get_mutant_name_str,
    get_involved_mutation,
    is_mutant_wt,
)
