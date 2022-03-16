from typing import Any

class MultiwfnConfig:
    def __init__(self, parent=None):
        self._parent = parent
    
    def required_executables(self):
        return [ self.EXE ]
    
    def required_env_vars(self):
        return [ self.DIR ]
    
    # -----------------------------
    # Cores for Multiwfn job (higher pirority)
    #
    # n_cores = n_cores
    # -----------------------------
    # Per core memory in MB for Multiwfn job (higher pirority)
    #
    # max_core = max_core
    # -----------------------------
    # Executable Multiwfn command for current environment
    #
    EXE = "Multiwfn"
    # -----------------------------
    # Path of Multiwfn folder
    #
    DIR = "$Multiwfnpath"
    
    
    def __getitem__(self, key : str ) -> Any:
        return getattr(self,key)
     
    def __setitem__(self, key : str, value: Any ) -> None:
        is_path = False
        if key in {"DIR", "EXE"}:
            is_path = True
        
        setattr(self,key,value)
        if is_path:
            self._parent.update_paths()
    
