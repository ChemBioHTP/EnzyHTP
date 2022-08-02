"""
quantum_mechanics module for enzy_htp. Allows for interfacing with quantum mechanics 
packages to analyze and manipulate structures. Follows pattern of having two companion classes
for each software package. First, a <Package>Interface() which uses the external software, and
second a <Package>Config() class which holds configuration settings for that Package. Each 
<Packge>_config submodule also contains a default_<Package>_config method which generates 
a default <Package>Config() describing the default behavior seen within enzy_htp.

Supported packages include:
	+ Gaussian 
	+ Multiwfn

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-11
"""
#from .gaussian_interface import GaussianInterface
#from .gaussian_config import GaussianConfig, default_gaussian_config
#from .multiwfn_interface import MultiwfnInterface
#from .multiwfn_config import MultiwfnConfig, default_multiwfn_config

from .interface import Interface
