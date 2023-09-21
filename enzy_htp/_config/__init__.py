"""This is an implementation module for accessing configuration setings for external 
software accesion in EnzyHTP. Works as a companion to the enzy_htp._interface.Interface
class, supplying it with settings. SHOULD NOT be hanlded directly by users. There is 
instead a singleton attribute enzy_htp.config with an instance. 

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-12
"""
from typing import List

from .config import Config

config = Config()
"""
Singleton object for accessing all configurations.
"""

import os
from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER

config_file: str = os.path.expandvars('$HOME/.eh_config')

if fs.has_content(config_file):
    _LOGGER.info(f"Found config file: {config_file}...")
    lines: List[str] = fs.lines_from_file(config_file)
    counter: int = 0
    for ll in lines:
        tks = ll.split('#')
        tk = tks[0].strip()
        if not tk:
            continue
        exec(tk)
        counter += 1

    _LOGGER.info(f"Updated {counter} config settings!")
