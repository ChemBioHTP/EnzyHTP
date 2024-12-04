"""the submodule for operation of Structure objects. This module should just rely on Structure and below
so that it can be used by anywhere that uses Structure.
It is seperated from the Structure class just for maintenance need. (i.e.: seperate the DS code and the operation code)

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2022-09-19
"""
from .general import (
    remove_solvent,
    remove_counterions,
    remove_hydrogens,
    remove_empty_chain,
    remove_non_peptide,
    update_residues,
    align_atom_order_in_each_residue,
)
from .residue import (
    deprotonate_residue,
    get_default_deproton_info,
    remove_side_chain_mutating_atom,
    check_res_topology_error,
    closest_n_residues
)
