"""This is an implementation module for accessing configuration setings for external 
software accesion in EnzyHTP. Works as a companion to the enzy_htp._interface.Interface
class, supplying it with settings. SHOULD NOT be hanlded directly by users. There is 
instead a singleton attribute enzy_htp.config with an instance.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-12
"""

from .config import Config

config = Config()
"""
Singleton object for accessing all configurations.
"""
