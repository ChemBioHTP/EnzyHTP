"""Module for general interface of preparation module. Aim for simply the use
of this module for light users.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-10-21
"""

from typing import List
from enzy_htp.structure import Structure


def detect_prepare_problem(stru: Structure, **kwargs) -> List[str]:
    """
    place holder for the integrate automatic detection of preparation
    problems.
    """


def prepare_stru(stru: Structure, **kwargs) -> Structure:
    """
    place holder for the integrate automatic for preparing the structure with
    detected problems.
    """
