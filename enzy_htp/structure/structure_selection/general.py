"""The module define functions that select a portion of the structure with a pymol-like pattern

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-15
"""
from enzy_htp.structure import Structure

def select_stru(stru: Structure, pattern: str) -> StruSelection:
    """use a pymol-like pattern to select part of the {stru}
    Args:
        stru: the Structure object of reference
        pattern: a pymol-like syntax to select residue positions
    Returns:
        a StruSelection object containing
        **reference** of the selection part of the {stru}"""
    pass

