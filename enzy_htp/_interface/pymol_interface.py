"""Defines an PyMolInterface class that serves as a bridge for enzy_htp to utilize PyMol package.
What does it need besides just functions realizing pymol functions?
- configuration (why cant it just be reading constant from another module?)
- environmental check (this should be enforced only when the interface is used)
The pymol2 python package is a wrapper around Cpp code. So it is not possible to
directly convert data structures but can only use pymol2 provided cmd.xxx APIs.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-02-15
"""

from typing import Tuple, List

from enzy_htp import config as project_config
from enzy_htp.core import env_manager as em
from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.pymol_config import PyMolConfig, default_pymol_config
from enzy_htp.structure import Structure, PDBParser

# QZ:Env check for package type interface?
try:
    import pymol2
except ImportError:
    _LOGGER.error("pymol package is missing. Install with `conda install -c conda-forge -c schrodinger pymol-bundle`")

class PyMolInterface:
    """Class that provides a direct inteface for enzy_htp to utilize PyMol package. 
    Atributes:
        config_ : the PyMolConfig() class which provides settings of PyMol 
    """

    def __init__(self, config: PyMolConfig = None) -> None:
        """Simplistic constructor that optionally takes an PyMolConfig object as its only
        argument.
        """
        self.config_ = config
        if not self.config_:
            self.config_ = default_pymol_config()

    def new_pymol_session(self) -> pymol2.PyMOL:
        """create a new pymol session, start it, and set feedback level
        TODO(qz): is making a pymol_session class a good idea?"""
        result = pymol2.PyMOL()
        result.start()
        result.cmd.feedback(*self.config.DEFAULT_OUTPUT_LV)
        return result

    # == intra-session modular functions == (requires a session and wont close it)
    def load_enzy_htp_stru(self, stru: Structure, pymol_session: pymol2.PyMOL) -> Tuple[str, pymol2.PyMOL]:
        """convert enzy_htp.Structure into a pymol object in a
        pymol2.PyMOL() session. 
        Using cmd.load and mediate by PDB format
        
        Returns:
            (pymol_obj_name, pymol_session) since the name is only valid in the session
            
        Note: pymol wont reset residue idx or chain names but will reset atom index from 1"""

        # create temp PDB
        pdb_str = PDBParser().get_file_str(stru, if_renumber=False)
        temp_dir = project_config["system.SCRATCH_DIR"]
        temp_pdb_path = fs.get_valid_temp_name(f"{temp_dir}/temp_pymol_interface.pdb")
        fs.safe_mkdir(temp_dir)
        with open(temp_pdb_path, "w") as f:
            f.write(pdb_str)
        # create pymol obj
        pymol_obj_name = pymol_session.cmd.get_unused_name("enzy_htp_stru")
        pymol_session.cmd.load(temp_pdb_path, pymol_obj_name)

        # clean up temp PDB
        if _LOGGER.level > 10:  # not DEBUG or below
            fs.safe_rm(temp_pdb_path)
            fs.safe_rmdir(temp_dir, empty_only=True)

        return (pymol_obj_name, pymol_session)

    def select_pymol_obj(self, pattern: str, pymol_obj_name: str, pymol_session: pymol2.PyMOL) -> List[int]:
        """an internal function return atom indexes of a pymol selection of a pymol obj 
        in the pymol session
        NOTE: in distance related selection scheme, there should be only on object in the session

        Returns:
            a list of atom indexes 
            if the object is created from enzy_htp.Structure(), these indexex will be
            correponding to Atom().idx in Structure()"""
        result = []
        # san check
        if len(pymol_session.cmd.get_object_list()) != 1:
            _LOGGER.warning(
                "more than 1 object is found in current pymol session when doing select. May give absurd results if distance based selection is used. (It will use atoms from other objs)")
        sele_pattern = f"{pymol_obj_name} & ({pattern})"
        sele_name = pymol_session.cmd.get_unused_name("enzy_htp_stru_sele")
        pymol_session.cmd.select(sele_name, sele_pattern)
        pymol_session.cmd.iterate(sele_name, "result.append(ID)", space=locals()) # use ID instead of index here

        return sorted(result)

    # == inter-session modular functions == (do not requires a session, will start and close one)
    # pass

    # TODO: go to a parent class
    @property
    def config(self) -> PyMolConfig:
        """Getter for the PyMolConfig() instance belonging to the class."""
        return self.config_

    def display_config(self) -> None:
        """Prints all settings for the object's PyMolConfig() inteface to stdout using PyMolConfig.display()."""
        self.config_.display()

    #TODO support pymol2 for these
    # def get_charge(self, fname: str, sele: str = '(all)') -> int:
    #     """Method that gets the formal charge for the specified sele in the 
    #     specified file. File must be a supported file type as listed in 
    #     PyMOLConfig. Checks if file exists and is supported type. If either 
    #     are not true then an error is logged and the program exits. DOES
    #     NOT check if the sele is valid. 
    #     Args:
    #         fname: The str() path to the file in question.
    #         sele: The str() specifying the atom selection in PyMOL synatx. Default is '(all)'.
    #     Returns:
    #         The formal charge of the sele as an int().
    #     """
    #     if not Path(fname).exists():
    #         _LOGGER.error(f"The supplied file '{fname}' does not exist. Exiting...")
    #         exit(1)

    #     if not self.supported_file_type(fname):
    #         _LOGGER.error(
    #             f"The supplied file '{fname}' has an unsupported extension. Exiting...")
    #         exit(1)

    #     _eh_local: Dict[str, Any] = {'fc': []}

    #     self.cmd.delete('all')
    #     self.cmd.load(fname)
    #     self.cmd.iterate(sele, 'fc.append(formal_charge)', space=_eh_local)
    #     self.cmd.delete('all')

    #     return sum(_eh_local['fc'])

    # def get_sequence(self, fname: str, sele: str = '(all)') -> str:
    #     """Method that gets the sequence for the specified sele in the 
    #     specified file. File must be a supported file type as listed in 
    #     PyMOLConfig. Checks if file exists and is supported type. If either 
    #     are not true then an error is logged and the program exits. DOES
    #     NOT check if the sele is valid. Note that non-canonical residues will
    #     be represented by a '?' in the sequence.
    #     Args:
    #         fname: The str() path to the file in question.
    #         sele: The str() specifying the atom selection in PyMOL synatx. Default is '(all)'.
    #     Returns:
    #         A str() with the sequence of the sele as it would appear in a .fasta file..
    #     """
    #     if not Path(fname).exists():
    #         _LOGGER.error(f"The supplied file '{fname}' does not exist. Exiting...")
    #         exit(1)

    #     if not self.supported_file_type(fname):
    #         _LOGGER.error(
    #             f"The supplied file '{fname}' has an unsupported extension. Exiting...")
    #         exit(1)

    #     self.cmd.delete('all')
    #     self.cmd.load(fname)
    #     lines: List[str] = self.cmd.get_fastastr(sele).splitlines()
    #     self.cmd.delete('all')

    #     lines = list(filter(lambda ll: ll[0] != '>', lines))
    #     return ''.join(lines)

class OpenPyMolSession:
    """a context manager that open a pymol session once enter and close once exit"""
    def __init__(self, pymol_interface: PyMolInterface) -> None:
        self.interface = pymol_interface

    def __enter__(self):
        """open a pymol session once enter"""
        self.session = self.interface.new_pymol_session()
        return self.session

    def __exit__(self, exc_type, exc_val, exc_tb):
        """open a pymol session once enter"""
        self.session.stop()
