from typing import Any
from plum import dispatch

#TODO documentation
class AmberConfig:
    def __init__(self, parent=None ):
        self._parent = parent

    def required_executables(self):
        return [ self.CPU_ENGINE, self.GPU_ENGINE, 'tleap', 'ampdb' ]

    def required_env_vars(self):
        return [ self.HOME ]

    # -----------------------------
    # Amber Home dir
    #
    HOME = "$AMBERHOME"
    # -----------------------------
    # User defined Amber CPU path (default: $AMBERHOME/bin/sander.MPI)
    #
    CPU_ENGINE = "$AMBERHOME/bin/sander.MPI"
    # -----------------------------
    # User defined Amber GPU path (default: $AMBERHOME/bin/pmemd.cuda)
    #
    GPU_ENGINE = "$AMBERHOME/bin/pmemd.cuda"
    # -----------------------------
    # Default configuration for MD
    #
    # -----------------------------
    # Water box type and size
    # availible type: box/oct
    BOX_TYPE = "oct"
    BOX_SIZE = "10"
    # -----------------------------
    #            min
    # Minimize
    #  &cntrl
    #   imin  = 1,  ntx   = 1,  irest = 0,
    #   ntc   = 2,    ntf = 2,
    #   cut   = 10.0,
    #   maxcyc= 20000, ncyc  = {0.5maxcyc},
    #   ntpr  = {0.01maxcyc},   ntwx  = 0,
    #   ntr   = 1,	restraint_wt = 2.0, restraintmask = '@C,CA,N',
    #  /
    CONF_MIN = {
        "ntc": "2", "ntf": "2",
        "cut": "10.0",
        "maxcyc": 20000,
        "ncyc": "0.5maxcyc",
        "ntpr": "0.01maxcyc",
        "ntr": "1",
        "restraintmask": "'@C,CA,N'",
        "restraint_wt": "2.0",
    }
    # -----------------------------
    #            heat
    # Heat
    #  &cntrl
    #   imin  = 0,  ntx = 1, irest = 0,
    #   ntc   = 2, ntf = 2,
    #   cut   = 10.0,
    #   nstlim= 20000, dt= 0.002,
    #   tempi = 0.0,  temp0=300.0,
    #   ntpr  = {0.01nstlim},  ntwx={nstlim},
    #   ntt   = 3, gamma_ln = 5.0,
    #   ntb   = 1,  ntp = 0,
    #   iwrap = 1,
    #   nmropt= 1,
    #   ig    = -1,
    #   ntr   = 1,	restraint_wt = 2.0, restraintmask = '@C,CA,N',
    #  /
    #  &wt
    #   type  = 'TEMP0',
    #   istep1= 0, istep2={0.9nstlim},
    #   value1= 0.0, value2=300.0,
    #  /
    #  &wt
    #   type  = 'TEMP0',
    #   istep1= {A_istep2+1}, istep2={nstlim},
    #   value1= 300.0, value2=300.0,
    #  /
    #  &wt
    #   type  = 'END',
    #  /
    CONF_HEAT = {
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
    # -----------------------------
    #            equi
    # Equilibration: constant pressure
    #  &cntrl
    #   imin  = 0,  ntx = 5,  irest = 1,
    #   ntf   = 2,  ntc = 2,
    #   nstlim= 500000, dt= 0.002,
    #   cut   = 10.0,
    #   temp0 = 300.0,
    #   ntpr  = {0.002nstlim}, ntwx = 5000,
    #   ntt   = 3, gamma_ln = 5.0,
    #   ntb   = 2,  ntp = 1,
    #   iwrap = 1,
    #   ig    = -1,
    #   ntr   = 1,	restraint_wt = 2.0,
    #   restraintmask = '@C,CA,N',
    #  /
    CONF_EQUI = {
        "ntx": "5",
        "irest": "1",
        "ntc": "2",
        "ntf": "2",
        "cut": "10.0",
        "nstlim": 500000,
        "dt": "0.002",
        "temp0": "300.0",
        "ntpr": "0.002nstlim",
        "ntwx":
            "5000",  # default 10ps (TODO support different power numbers)
        "ntt": "3",
        "gamma_ln": "5.0",
        "iwarp": "1",
        "ntr": "1",
        "restraintmask": "'@C,CA,N'",
        "restraint_wt": "2.0",  # the later two are only used when ntr = 1
    }
    # -----------------------------
    #            prod
    # Production: constant pressure
    #  &cntrl
    #   imin  = 0, ntx = 1, irest = 0,
    #   ntf   = 2,  ntc = 2,
    #   nstlim= 50000000, dt= 0.002,
    #   cut   = 10.0,
    #   temp0 = 300.0,
    #   ntpr  = 50000, ntwx = 5000,
    #   ntt   = 3, gamma_ln = 5.0,
    #   ntb   = 2,  ntp = 1,
    #   iwrap = 1,
    #   ig    = -1,
    #  /
    CONF_PROD = {
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

    def __getitem__(self, key : str ) -> Any:
        return getattr(self,key)
     
    def __setitem__(self, key : str, value: Any ) -> None:
        is_path = False
        if key in {"CPU_ENGINE", "GPU_ENGINE", "HOME"}:
            is_path = True
        
        setattr(self,key,value)
        if is_path:
            self._parent.update_paths()

  
