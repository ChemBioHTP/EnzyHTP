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

from .load_config import load_config

import os
from enzy_htp.core import file_system as fs


config_file: str = os.path.expandvars('$HOME/.eh_config')


if fs.has_content(config_file):
    config.load_config( config_file )

