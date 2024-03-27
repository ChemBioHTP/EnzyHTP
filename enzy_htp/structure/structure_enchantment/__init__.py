"""the submodule for initialization of atomic properties of Structure objects. (e.g.: charge, connectivity)
this module is importing _interface

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2024-03-26
"""
from .connectivity import (
    init_connectivity
)

from .charge import init_charge
