"""
TODO
"""

from enzy_htp import core
from .amber_interface import AmberInterface
from .gaussian_interface import GaussianInterface
from .multiwfn_interface import MultiwfnInterface

class Interface:


    def __init__(self, config):
        self.amber = AmberInterface( config._amber )
        self.gaussian = GaussianInterface( config._gaussian )
        self.multiwfn = MultiwfnInterface( config._multiwfn )

