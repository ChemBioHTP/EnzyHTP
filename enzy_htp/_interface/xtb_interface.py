"""Defines an XTBInterface class that serves as a bridge for enzy_htp to utilize the xtb
semi-empirical quantum mechanics (QM) package. Uses the XTBConfig class found in enzy_htp/_config/xtb_config.py.
Supported operations include:

    + single point energy calculations
    + geometry optimization 


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-08-17
"""
from pathlib import Path
from typing import Tuple, List

from .base_interface import BaseInterface

from enzy_htp import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.xtb_config import XTBConfig, default_xtb_config


class XTBInterface(BaseInterface):
    """Class that provdes a direct interface for enzy_htp to utilize the xtb semi empirical method. Supported methods
    include single energy calculations and geometry optimization
    
    Attributes:
        config_ : The XTBXonfig() class which provides settings for both running xtb and maintaing a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : a bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config: XTBConfig = None) -> None:
        """Simplicstic constructor that requires the parent interface as an argument and optionally takes an XTBConfig instance.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_xtb_config)

    def single_point(self, fname: str, charge: int = 0, spin: int = 1, n_iter: int = -1, n_proc: int = -1, work_dir:str=None) -> float:
        """Performs a single point energy calculation on the supplied file. Checks if the file exists and is one of the correct
        file formats. Errors if not. Returns a value in Hartrees. 

        Args:
            fname: Name of the input file to use.
            charge: Optional charge of the system as an int(). Defaults to 0.
            spin: Optional spin of the system as an int(). Defaults to 0.
            n_iter: Number of SCF iterations. Optional, if not specified uses value in XTBConfig.
            n_proc: Number of processes to use in calculation. Optional, if not specified uses value in XTBConfig.
            work_dir:TODO(CJ)

        Returns:
            Single point energy value in Hartrees.           
        """

        fs.check_file_exists(fname)
        self._check_valid_extension(fname)

        if n_iter == -1:
            n_iter = self.config_.N_ITER

        if n_proc == -1:
            n_proc = self.config_.N_PROC

        results = self.env_manager_.run_command(
            self.config_.XTB_EXE,
            ["--chrg", str(charge), "--iterations",
             str(n_iter), "--parallel",
             str(n_proc), "--norestart", "--sp", fname])

        self._remove_temp_files(str(Path(fname).parent))

        for ll in results.stdout.splitlines():
            if ll.find('TOTAL ENERGY') != -1:
                return float(ll.split()[3])

        _LOGGER.error(f"ERROR: Single point energy calculation for file '{fname}' did not contain energy value. Exiting...")
        exit( 1 )

    def geo_opt(self,
                fname: str,
                charge: int = 0,
                spin: int = 1,
                n_iter: int = -1,
                n_proc: int = -1,
                freeze_atoms:List[int]=None,
                constraints:List=None,
                work_dir:str=None
                ) -> Tuple[str, float]:
        """TODO(CJ)

        Args:

        Returns:
            
        """
        fs.check_file_exists(fname)
        self._check_valid_extension(fname)

        fs.check_file_exists(fname)
        self._check_valid_extension(fname)

        if n_iter == -1:
            n_iter = self.config_.N_ITER

        if n_proc == -1:
            n_proc = self.config_.N_PROC

        input_lines:List[str]=list()

        if freeze_atoms is not None:
            input_lines.extend([
                '$fix',
               f'   atoms: {",".join(map(str,freeze_atoms))}',
                '$end'
            ])

        if constraints is not None:
            input_lines.append('$constrain')
            for cst in constraints:
                input_lines.append(
                f'   {cst[0]}: {", ".join(map(lambda ff: str(int(ff)), cst[-1]))}, {cst[1]:.2f}'
                )
            input_lines.append('$end')

        cmd_args:List[str] = ["--chrg", str(charge), "--iterations", str(n_iter), "--parallel", str(n_proc), "--norestart" ]

        geom_input:str=None
        if input_lines:
            temp_path = Path(fname)
            geom_input = str(temp_path.parent / f"__geo_opt.inp")
            fs.write_lines(geom_input, input_lines)
            cmd_args.extend(["--input", geom_input ]) 

        cmd_args.extend([fname, "--opt" ])
        ext:str=fs.get_file_ext(fname)
        expected_output:str=f"xtbopt{ext}"

        results = self.env_manager_.run_command(
                self.config_.XTB_EXE,
                cmd_args)
        
        self._remove_temp_files( str(Path(fname).parent) )
        
        energy:float=None
        for ll in results.stdout.splitlines():
            if ll.find('TOTAL ENERGY') != -1:
                energy = float(ll.split()[3])

        self._remove_temp_files("./")
    
        return (expected_output, energy)


    def _remove_temp_files(self, work_dir: str, extra_files:List[str]=None) -> None:
        """Removes the expected temp files created in the working directory of the input file. Deletes the file safely and 
        silently.

        Args:
            work_dir: Name of the directory to search for the files in.
            extra_files: A List[str] containing extra files to delete in the function.

        Returns:
            Nothing.
        """

        for fname in "charges wbo xtbrestart xtbtopo.mol xtbopt.log xtbopt.pdb xtbopt.mol xtbopt.xyz .xtboptok".split():
            fs.safe_rm(f"{work_dir}/{fname}")

        if extra_files:
            for ef in extra_files:
                fs.safe_rm(ef)

    def _check_valid_extension(self, fname: str) -> None:
        """Does the supplied file have a supported file extension? Logs an error
        and exits if not."""
        ext: str = fs.get_file_ext(str(fname))
        if ext in self.config_.SUPPORTED_EXTENSIONS:
            return

        _LOGGER.error(
            f"The supplied file '{fname}' has an unsupported extension of '{ext}'. Supported include {', '.join(self.config_.SUPPORTED_EXTENSIONS)}. Exiting..."
        )
        exit(1)
