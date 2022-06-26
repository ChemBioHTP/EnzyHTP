"""
quantum_mechanics module for enzy_htp. Allows for interfacing with quantum mechanics 
packages to analyze and manipulate structures. Follows pattern of having two companion classes
for each software package. First, a <Package>Interface() which uses the external software, and
second a <Package>Config() class which holds configuration settings for that Package.

Supported packages include:
	+ Gaussian 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-11
"""
from .gaussian_interface import GaussianInterface
from .gaussian_config import GaussianConfig, default_gaussian_config
