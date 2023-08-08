"""The module define functions that select a portion of the structure with a pymol-like pattern

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-15
"""
from enzy_htp import interface, OpenPyMolSession
from enzy_htp.structure import Structure
from .stru_selection import StruSelection


def select_stru(stru: Structure, pattern: str) -> StruSelection:
    """use a pymol-like pattern to select part of the {stru}
    use pymol selection as engine. (it is better to have only one engine
    so this function is here instead of pymol interface)
    TODO(qz): in the future when it is necessary we should parse the selection
              by enzyhtp itself.
    Args:
        stru: the Structure object of reference
        pattern: a pymol-like syntax to select residue positions
    Returns:
        a StruSelection object containing
        **reference** of the selection part of the {stru}

    Details:
        1. convert Structure to a pymol object (in pymol interface)
        2. use sele in pymol get a list of atom idxes (in pymol interface)
        3. get these atoms and construct a StruSelection object from them"""
    pi = interface.pymol
    with OpenPyMolSession(pi) as pms:
        pymol_obj_name = pi.load_enzy_htp_stru(stru, session=pms)[0]
        atom_idx_list = pi.select_pymol_obj(pattern, pymol_obj_name, pms)
    # construct a selection object base on a list of atom indexes
    return StruSelection.from_atom_idx_list(stru, atom_idx_list)
