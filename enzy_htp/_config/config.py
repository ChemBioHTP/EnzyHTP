"""Sub module that lays out the Config() class which stores all configuration settings
for the EnzyHTP package. Config() is meant to be used as a singleton and typically as the
enzy_htp.config attribute. DO NOT create instances of this object on your own. It is created
only by the module. For each <Package> that EnzyHTP interfaces with, there should be a
<Package>Interface in enzy_htp/_interface and a <Package>Config in enzy/_config. At present,
there are configuration settings for the below packages by the given <Package>Config:

    + SystemConfig
    + ARMerConfig
    + AlphaFill, AlphaFillConfig 
    + AmberMD, AmberConfig
    + BCL, BCLConfig
    + Gaussian, GaussianConfig
    + MOE, MOEConfig
    + Mole2, Mole2Config
    + Multiwfn, MultiwfnConfig
    + PyMOL, PyMolConfig
    + Rosetta, RosettaConfig
    + RDKit, RDKitConfig
    + xtb, XTBConfig

In addition to specific packages, settings for the system are specified with SystemConfig

Note that for a given <Package>, there should exist a default_<package>_config function which
provides a default, deep-copied instance of the specific <Package>Confg.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-07-15
"""
from typing import Any, Dict, List
from .alphafill_config import AlphaFillConfig, default_alphafill_config
from .amber_config import AmberConfig, default_amber_config
from .bcl_config import BCLConfig, default_bcl_config
from .gaussian_config import GaussianConfig, default_gaussian_config
from .moe_config import MOEConfig, default_moe_config
from .mole2_config import Mole2Config, default_mole2_config
from .multiwfn_config import MultiwfnConfig, default_multiwfn_config
from .pymol_config import PyMolConfig, default_pymol_config
from .rdkit_config import RDKitConfig, default_rdkit_config
from .rosetta_config import RosettaConfig, default_rosetta_config
from .system_config import SystemConfig, default_system_config
from .armer_config import ARMerConfig, default_armer_config
from .xtb_config import XTBConfig, default_xtb_config

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs


class Config:
    """Class that holds all configuration settings for the different external softwares 
    that EnzyHTP relies on. Uses a parameter setting/getting syntax that uses [] operators
    as well as namespace dotting. Direct accession of the sub-configs is NOT recommended. 

    Attributes:
        _alphafill: Private instance of AlphaFillConfig() with default settings.
        _amber: Private instance of AmberConfig() with default settings.
        _bcl: Private instance of BCLConfig() with default settings.
        _gaussian: Private instance of GaussianConfig() with default settings.
        _moe: Private instance of MOEConfig() with default settings.
        _mole2: Private instance of Mole2Config() with default settings.
        _multiwfn: Private instance of MultiwfnConfig() with default settings.
        _pymol: Private instance of PyMolConfig() with default settings.
        _rdkit: Private instace of RDKitConfig() with default settings.
        _rosetta: Private instance of RosettaConfig() with default settings.
        _system: Private instance of SystemConfig() with default settings.
        _armer: Private instance of ARMerConfig() with default settings.
        _xtb: Private instance of XTBConfig() with default settings.
        _mapper: Helper dict() that maps between name of the package and its config. For internal use only.
    """

    def __init__(self):
        """Constructor that creates a <Package>Config instance for each <package> using default_<package>_config."""
        self._alphafill = default_alphafill_config()
        self._amber = default_amber_config()
        self._bcl = default_bcl_config()
        self._gaussian = default_gaussian_config()
        self._multiwfn = default_multiwfn_config()
        self._moe = default_moe_config()
        self._mole2 = default_mole2_config()
        self._pymol = default_pymol_config()
        self._rdkit = default_rdkit_config()
        self._rosetta = default_rosetta_config()
        self._system = default_system_config()
        self._armer = default_armer_config()
        self._xtb = default_xtb_config()
        self._mapper = {
            "alphafill": self._alphafill,
            "amber": self._amber,
            "bcl": self._bcl,
            "gaussian": self._gaussian,
            "moe": self._moe,
            "mole2": self._mole2,
            "multiwfin": self._multiwfn,
            "pymol": self._pymol,
            "rdkit": self._rdkit,
            "rosetta": self._rosetta,
            "system": self._system,
            "xtb": self._xtb,
            "armer": self._armer,
        }

    def __getitem__(self, key: str) -> Any:
        """Getter for the settings in the Config() object. Uses the grammar: "<package>.<setting>" 
        Note that the <packge> is specified in all lowercase letters.

        Args:
            key: A string with the grammar "<package>.<name>" that you will returning.

        Returns:
            Corresponding value, if it exists.

        """
        if key.count("."):
            app, settings = key.split(".", 1)
            ptr = self._mapper.get(app, None)
            if ptr is not None:
                return ptr[settings]

        _LOGGER.error(f"The supplied key {key} is invalid. Exiting...")
        exit(1)

    def __setitem__(self, key: str, value: Any) -> None:
        """Setter for the configuration values in the Config object. Uses the grammar
        config[key] = value, where key has the form "<package>.<name>" and value can be 
        any value. Note that the <packge> is specified in all lowercase letters.

        Args:
            key: A string with the grammar "<package>.<name>" for the value you will be updating.
            value: The value you are going to set the corresponding variable to. Can have any type.

        Returns:
            Nothing.

        """
        if key.count("."):
            app, settings = key.split(".", 1)
            ptr = self._mapper.get(app, None)
            if ptr is not None:
                ptr[settings] = value
                return

        _LOGGER.error(f"The supplied key {key} is invalid. Exiting...")
        exit(1)


    def load_config(self, fname: str) -> None:
        """Function that loads configuration settings from the supplied file. The supplied file
        can have comments. The function utilizes basic type deduction to interpret the supplied settings.

        Args:
            fname: Name of the file to parse a str(). 

        Returns:
            Nothing.
        """
        if not fs.has_content(fname):
            _LOGGER.error(f"The supplied file {fname} does not exist or is empty. Exiting...")
            exit( 1 )

        _LOGGER.info(f"Found config file: {fname}...")
        lines: List[str] = fs.lines_from_file(fname)
        
        counter: int = 0
        for ll in lines:
            tks = ll.split('#')
            tk = tks[0].strip()
            if not tk:
                continue
            
            key,value=tk.split('=',1)
            key,value=key.strip(),value.strip()
            if key[0:6] != 'config':
                _LOGGER.warning(f"Improper syntax found in line: {tk}. Continuing...")
                continue
            key = key[7:-1]
            
            key = key.replace("\'","")
            key = key.replace("\"","")

            if value in ["True", "False"]:
                value = bool(value)                
            elif value == "None":
                value = None
            elif value[0].isnumeric():
                if value.find('.') != -1:
                    value = float(value)
                else:
                    value = int(value)
            
            elif value[0] == value[-1] and value[0] in ["\"", "\'"]:
                value = value.replace("\'","")
                value = value.replace("\"","")
            else:
                _LOGGER.warning(f"Improper syntax found in value {value}. Continuing...")
                continue
            
            self[key] = value
            counter += 1

        _LOGGER.info(f"Updated {counter} config settings!")

