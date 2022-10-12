"""This is an implementation module for interfacing with external software packages
in EnzyHTP. Works as a companion to the enzy_htp._config.Config class, using the former
as a source of settings. SHOULD NOT be hanlded directly by users. There is 
instead a singleton attribute enzy_htp.config with an instance.
Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-12
"""

from .interface import Interface