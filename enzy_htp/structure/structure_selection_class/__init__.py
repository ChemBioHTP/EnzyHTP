"""this module is created to avoid loop importing issue. The fact that in
importing layers, structure_selection APIs uses PyMol from interface while
the datastructure StruSelection is used by Amber, Gaussian, from interface.

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-12-26
"""

from .stru_selection import StruSelection