"""Defines AmberConfig() which holds configuration settings for enzy_htp to interface with the 
Amber software package. In addition to holding Amber settings, AmberConfig() creates the input
files for minimization, heating, constant pressure production, and constant pressure
equilibration. File also contains default_amber_config() which creates a default version
of the AmberConfig() object.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-02
"""
from copy import deepcopy
from typing import Any, List, Dict
from plum import dispatch


class AmberConfig:
    """Class that holds default values for running Amber within enzy_htp and also creates
    input files for minimzation, heating, constant pressure production, and constant
    pressure equilibration.

    Attributes:
            parent_ : Points to parent config object. Optional and defaults to None.
            HOME : str() corresponding to Amber home directory on the system.
            CPU_ENGINE : str() corresponding to Amber cpu sander.
            GPU_ENGINE : str() corresponding to Amber gpu sander.
            BOX_TYPE : str() corresponding to type of water box.
            BOX_SIZE : TODO
            CONF_MIN : dict() holding settings for Amber minimization.
            CONF_HEAT : dict() holding settings for Amber heating.
            CONF_EQUI : dict() holding settings for Amber constant pressure equilibration run.
            CONF_PROD : dict() holding settings for Amber constant pressure production run.
    """

    HOME: str = "$AMBERHOME"
    """Environment variable for Amber's HOME Directory"""

    CPU_ENGINE: str = "$AMBERHOME/bin/sander.MPI"
    """Environment variable for Amber's cpu engine"""

    GPU_ENGINE: str = "$AMBERHOME/bin/pmemd.cuda"
    """Environment variable for Amber's gpu engine"""

    BOX_TYPE: str = "oct"
    """Water Box type. Allowed values are 'box' and 'oct'"""
    # TODO(CJ): Put in checker for this
    BOX_SIZE = "10"
    """Water Box size."""
    # TODO(CJ): Should this be a str(), int() or float()?

    CONF_MIN: Dict = {
        "ntc": "2",
        "ntf": "2",
        "cut": "10.0",
        "maxcyc": 20000,
        "ncyc": "0.5maxcyc",
        "ntpr": "0.01maxcyc",
        "ntr": "1",
        "restraintmask": "'@C,CA,N'",
        "restraint_wt": "2.0",
    }
    """dict() holding the settings for an Amber minimization."""

    CONF_HEAT: Dict = {
        "ntc": "2",
        "ntf": "2",
        "cut": "10.0",
        "nstlim": 20000,
        "dt": "0.002",
        "tempi": "0.0",
        "temp0": "300.0",
        "ntpr": "0.01nstlim",
        "ntwx": "nstlim",
        "ntt": "3",
        "gamma_ln": "5.0",
        "iwarp": "1",
        "ntr": "1",
        "restraintmask": "'@C,CA,N'",
        "restraint_wt": "2.0",
        "A_istep2": "0.9nstlim",
        "B_istep1": "A_istep2+1",
    }
    """dict() holding the settings for an Amber heating."""

    CONF_EQUI: Dict = {
        "ntx": "5",
        "irest": "1",
        "ntc": "2",
        "ntf": "2",
        "cut": "10.0",
        "nstlim": 500000,
        "dt": "0.002",
        "temp0": "300.0",
        "ntpr": "0.002nstlim",
        "ntwx": "5000",  # default 10ps (TODO support different power numbers)
        "ntt": "3",
        "gamma_ln": "5.0",
        "iwarp": "1",
        "ntr": "1",
        "restraintmask": "'@C,CA,N'",
        "restraint_wt": "2.0",  # the later two are only used when ntr = 1
    }
    """dict() holding the settings for an Amber constant pressure equilibration run."""

    CONF_PROD: Dict = {
        "ntx": "5",
        "irest": "1",
        "ntc": "2",
        "ntf": "2",
        "cut": "10.0",
        "nstlim": 50000000,
        "dt": "0.002",
        "temp0": "300.0",
        "ntpr": "0.001nstlim",
        "ntwx": "5000",  # default 10ps
        "ntt": "3",
        "gamma_ln": "5.0",
        "iwarp": "1",
        "ntr": "0",
        "restraintmask": None,
        "restraint_wt": "2.0",  # the later two are only used when ntr = 1
    }
    """dict() holding the settings for an Amber constant pressure production run."""

    def __init__(self, parent=None):
        """Trivial constructor that optionally sets parent_ dependency. parent_ is None by default."""
        self.parent_ = parent

    def valid_box_type(self) -> bool:
        """Checks if the BOX_TYPE attribute is an acceptable value. Current allowed values are "box" and "oct"."""	 
        return self.BOX_TYPE in set("box oct".split())

    def required_executables(self) -> List[str]:
        """A hardcoded list of required executables for Amber."""
        return [self.CPU_ENGINE, self.GPU_ENGINE, "tleap", "ampdb"]

    def required_env_vars(self) -> List[str]:
        """A hardcoded list of required enviornment variables for Amber."""
        return [self.HOME]

    def __getitem__(self, key: str) -> Any:
        """Getter that enables [] accession of AmberConfig() attributes."""
        return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        """Setter that enables [] accession of AmberConfig() attributes with value validation."""
        print(key)
        is_path = False
        if key in {"CPU_ENGINE", "GPU_ENGINE", "HOME"}:
            is_path = True
        
        setattr(self, key, value)
        if not self.valid_box_type():
             raise TypeError() #TODO(CJ): make a custom error for this part

        if is_path and self.parent_:
            self.parent_.update_paths()

    def get_engine(self, mode: str) -> str:
        """Getter that returns the path to either the CPU or GPU engine configured for Amber.
		
		Args:
			mode: The mode that amber will be run in. Allowed values are "CPU" and "GPU".

		Returns:
			Path to the specified engine.

		Raises:
			TypeError if an invalid engine type is supplied.	
		"""
        if mode == "CPU":
            return self.CPU_ENGINE
        elif mode == "GPU":
            return self.GPU_ENGINE
        else:
			#TODO(CJ): add a custom error for this part
            raise TypeError()


def default_amber_config() -> AmberConfig:
    """Creates a deep-copied default version of the AmberConfig() class."""
    return deepcopy(AmberConfig())
