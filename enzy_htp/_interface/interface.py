"""Submodule describing an interface class which enables interfacing into all supported applications within EnzyHTP. This
model sees all interactions between EnzyHTP and a given application via a data attribute which is an instance of a 
specific application interface. At present, the below packages can be interacted with via the corresponding Interface
class. Packages:

    + AmberMD, AmberInterface
    + Gaussian, GaussianInterface
    + Multiwfn, MultiwfnInterface


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-20
"""

from enzy_htp import core
from .amber_interface import AmberInterface
from .gaussian_interface import GaussianInterface
from .multiwfn_interface import MultiwfnInterface

from enzy_htp._config import Config

class Interface:
    """Wrapper class that houses access to individual <Package>Interface classes that exist in EnzyHTP.
    Each <Package>Interface is available as the attribute Interface.<package> (all lower case). Each instance
    needs an EnzyHTP.Config instance to hold all data.

    Attributes:
        amber: Corresponds to instance of AmberInterface().
        gaussian: Corresponds to instnce of GaussianInterface().
        multiwfn: Corresponds to instance of MultiwfnInterface().
    """

    def __init__(self, config : Config):
        self.amber = AmberInterface(config._amber)
        self.gaussian = GaussianInterface(config._gaussian)
        self.multiwfn = MultiwfnInterface(config._multiwfn)
