"""


Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-07-15
"""
from typing import Any, Dict

from .amber_config import AmberConfig, default_amber_config
from .gaussian_config import GaussianConfig, default_gaussian_config
from .multiwfn_config import MultiwfnConfig, default_multiwfn_config
from .system_config import SystemConfig, default_system_config

from enzy_htp.core import _LOGGER


class Config:
    """Class that holds all configuration settings for the different external softwares 
    that EnzyHTP relies on. Uses a parameter setting/getting 


    Attributes:
        _amber: 
        _gaussian:
        _multiwfn:
        _system
    """
    def __init__(self):
        self._amber = default_amber_config()
        self._gaussian = default_gaussian_config()
        self._multiwfn = default_multiwfn_config()
        self._system  = default_system_config()

    def __getitem__(self, key: str ) -> Any:
        
        if key.count('.'):
            app, settings = key.split('.',1) 
            ptr = None
            if app == 'amber':
                ptr = self._amber
            elif app == 'gaussian':
                ptr = self._gaussian
            elif app == 'multiwfn':
                ptr = self._multiwfn
            elif app == 'system':
                ptr = self._system
            else:
                raise TypeError()
            return ptr[settings]
        else:
            raise TypeError()

    def __setitem__(self, key:str, value: Any) -> None:
        if key.count('.'):
            app, settings = key.split('.',1) 
            ptr = None
            if app == 'amber':
                ptr = self._amber
            elif app == 'gaussian':
                ptr = self._gaussian
            elif app == 'multiwfn':
                ptr = self._multiwfn
            elif app == 'system':
                ptr = self._system
            else:
                raise TypeError()
            ptr[settings] = value
        else:
            raise TypeError()

