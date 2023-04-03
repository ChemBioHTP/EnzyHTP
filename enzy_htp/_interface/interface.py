"""Submodule describing an interface class which enables interfacing into all supported applications within EnzyHTP. This
model sees all interactions between EnzyHTP and a given application via a data attribute which is an instance of a 
specific application interface. At present, the below packages can be interacted with via the corresponding Interface
class. Packages:
    + AmberMD, AmberInterface
    + BCL, BCLInterface
    + Gaussian, GaussianInterface
    + MOE, MOEInterface
    + Multiwfn, MultiwfnInterface
    + PyMOL, PyMOLInterface
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-07-20
"""

from enzy_htp import core
from .amber_interface import AmberInterface
from .bcl_interface import BCLInterface
from .gaussian_interface import GaussianInterface
from .moe_interface import MOEInterface
from .multiwfn_interface import MultiwfnInterface
from .pymol_interface import PyMOLInterface
from .rosetta_interface import RosettaInterface

from enzy_htp._config import Config


class Interface:
    """Wrapper class that houses access to individual <Package>Interface classes that exist in EnzyHTP.
    Each <Package>Interface is available as the attribute Interface.<package> (all lower case). Each instance
    needs an EnzyHTP.Config instance to hold all data.
    Attributes:
        amber: Corresponds to instance of AmberInterface().
        bcl: Corresponds to instance of BCLInterface().
        gaussian: Corresponds to instnce of GaussianInterface().
        moe: Corresponds to instance of MOEInterface().
        multiwfn: Corresponds to instance of MultiwfnInterface().
        pymol: Corresponds to instance of PyMOLInteface().
        rosetta: Corresponds to instance of RosettaInterface().
    """

    def __init__(self, config: Config):
        self.amber = AmberInterface(config._amber)
        self.bcl = BCLInterface(config._bcl)
        self.gaussian = GaussianInterface(config._gaussian)
        self.moe = MOEInterface(config._moe)
        self.multiwfn = MultiwfnInterface(config._multiwfn)
        self.pymol = PyMOLInterface(config._pymol)
        self.rosetta = RosettaInterface(config._rosetta)
