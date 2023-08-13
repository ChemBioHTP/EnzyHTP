"""the submodule for operation of Structure objects

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2022-09-19
"""
from .general import (
    remove_solvent,
    remove_empty_chain,
    remove_non_peptide,
    update_residues,
    order_atoms_to_stru,
)
from .residue import (
    deprotonate_residue,
    get_default_deproton_info,
    remove_side_chain_mutating_atom,
    check_res_topology_error,
)
