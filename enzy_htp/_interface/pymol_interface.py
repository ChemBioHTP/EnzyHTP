"""Defines an PyMolInterface class that serves as a bridge for enzy_htp to utilize PyMol package.
What does it need besides just functions realizing pymol functions?
- configuration (why cant it just be reading constant from another module?)
- environmental check (this should be enforced only when the interface is used)
The pymol2 python package is a wrapper around Cpp code. So it is not possible to
directly convert data structures but can only use pymol2 provided cmd.xxx APIs.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-02-15
"""
from pathlib import Path
from typing import List, Dict, Any, Union, Tuple

import pymol2
import pandas as pd

#TODO(CJ): add something to remove "PyMOL not running. Entering library mode (experimental" message on pymol running
from enzy_htp import config as eh_config
from enzy_htp.core import env_manager as em
from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER, check_var_type
from enzy_htp._config.pymol_config import PyMolConfig, default_pymol_config
from enzy_htp.structure import Structure, PDBParser
from enzy_htp.chemical.residue import THREE_LETTER_AA_MAPPER

from .base_interface import BaseInterface


class PyMolInterface(BaseInterface):
    """Class that provides a direct inteface for enzy_htp to utilize PyMol package. 
    
    Atributes:
        config_ : The PyMOLConfig() class which provides settings for both running PyMOL and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config: PyMolConfig = None) -> None:
        """Simplistic constructor that optionally takes an PyMOLConfig object as its only argument.
        Calls parent class.
        """
        super().__init__(parent, config, default_pymol_config)
        temp_session: pymol2.PyMOL = self.new_session()
        self.available_cmds_: List[str] = temp_session.cmd.__dict__['kw_list']

    def new_session(self) -> pymol2.PyMOL:
        """create a new pymol session, start it, and set feedback level. Canonical means to get a session in enzy_htp."""
        result = pymol2.PyMOL()
        result.start()
        result.cmd.feedback(*self.config.DEFAULT_OUTPUT_LV)
        return result

    def check_pymol2_installed(self) -> None:
        """Method that checks if """
        if "pymol2" in self.missing_py_modules():
            _LOGGER.errorf("pymol2 is NOT installed. Use 'conda install -c conda-forge -y -q pymol-open-source'. Exiting...")
            exit(1)

    def convert(self,
                session: pymol2.PyMOL,
                file_1: str,
                file_2: str = None,
                new_ext: str = None,
                split_states: bool = False) -> Union[str, List[str]]:
        """Method that converts a supplied file to a different format. Requires a pymol2 session from PyMolInterface.new_session().
        Either a new filename or new file  extension an be supplied. If neither or both are supplied, then the function will exist. 
        Note that the function does not check for valid file types and will catch any errors that are thrown if an invalid file_1 or
        output file combination is supplied. Returns the new outfile.

        Args:
            session : A pymol2.PyMOL() session to use
            file_1 : The name of the original file as a str().
            file_2 : The name of the output file as a str(). Optional.
            new_ext : The new extension to use. Optional.

        Returns:
            The name of the new file as a str().
        """
        self.check_pymol2_installed()
        #TODO(CJ): update split states. update return policy

        fs.check_file_exists(file_1)

        if (not file_2 and not new_ext):
            _LOGGER.error("Either a new file name or target extension must be supplied. Both are empty. Exiting...")
            exit(1)

        if file_2 and new_ext:
            _LOGGER.error("Either a new file name or target extension must be supplied. Both were supplied. Exiting...")
            exit(1)

        if not file_2:
            fpath = Path(file_1)
            file_2 = fpath.with_suffix(new_ext)

        result: str = str(file_2)

        try:
            session.cmd.delete('all')
            session.cmd.load(file_1)

            if split_states:
                session.cmd.split_states('all')
                result = Path(result)
                parent_dir = result.parent
                result_stem = result.stem
                result = []

                for oidx, oname in enumerate(session.cmd.get_object_list()):

                    if session.cmd.count_states(oname) != 1:
                        session.cmd.delete(oname)
                        continue

                    fname = str(parent_dir / f"{result_stem}_{oidx}{new_ext}")
                    session.cmd.save(fname, oname)
                    result.append(fname)
            else:
                session.cmd.save(result)
                session.cmd.delete('all')

        except Exception as exe:
            _LOGGER.error(f"Could not convert '{file_1}' to '{file_2}'. Encountered error: '{exe}'. Exiting...")
            exit(1)

        return result

    def supported_file_type(self, fname: str) -> bool:
        """Convenience function that checks if a listed filetype is supported by 
        the PyMOLInterface. Currently supported filetypes include:
            + .pdb
            + .mol2
            + .cif

        Args:
            fname: The str() path to the file in question.

        Returns:
            Whether the file type is supported.
        """
        #TODO(CJ): need to check which file formats are actually supported
        extension: str = Path(fname).suffix
        return extension in self.config().IO_EXTENSIONS

    def get_charge(self, fname: str, sele: str = '(all)', session: pymol2.PyMOL = None) -> int:
        """Method that gets the formal charge for the specified sele in the 
        specified file. File must be a supported file type as listed in 
        PyMOLConfig. Checks if file exists and is supported type. If either 
        are not true then an error is logged and the program exits. DOES
        NOT check if the sele is valid. 

        Args:
            fname: The str() path to the file in question.
            sele: The str() specifying the atom selection in PyMOL synatx. Default is '(all)'.
            session : A pymol2.PyMOL() session to use. By default uses an internal session instance that resets each function call.
        """
        pass

    # == intra-session modular functions == (requires a session and wont close it)
    def load_enzy_htp_stru(self, session: pymol2.PyMOL, stru: Structure) -> Tuple[str, pymol2.PyMOL]:
        """convert enzy_htp.Structure into a pymol object in a
        pymol2.PyMOL() session. 
        Using cmd.load and mediate by PDB format
        
        Returns:
            (pymol_obj_name, session) since the name is only valid in the session
            
        Note: pymol wont reset residue idx or chain names but will reset atom index from 1"""

        # create temp PDB
        self.check_pymol2_installed()
        pdb_str = PDBParser().get_file_str(stru, if_renumber=False)
        temp_dir = eh_config["system.SCRATCH_DIR"]
        temp_pdb_path = fs.get_valid_temp_name(f"{temp_dir}/temp_pymol_interface.pdb")
        fs.safe_mkdir(temp_dir)
        with open(temp_pdb_path, "w") as f:
            f.write(pdb_str)
        # create pymol obj
        pymol_obj_name = session.cmd.get_unused_name("enzy_htp_stru")
        session.cmd.load(temp_pdb_path, pymol_obj_name)

        # clean up temp PDB
        fs.clean_temp_file_n_dir([temp_dir, temp_pdb_path])

        return (pymol_obj_name, session)

    def select_pymol_obj(self, session: pymol2.PyMOL, pattern: str, pymol_obj_name: str) -> List[int]:
        """an internal function return atom indexes of a pymol selection of a pymol obj 
        in the pymol session
        NOTE: in distance related selection scheme, there should be only on object in the session

        Returns:
            a list of atom indexes 
            if the object is created from enzy_htp.Structure(), these indexex will be
            correponding to Atom().idx in Structure()"""
        self.check_pymol2_installed()
        result = []
        # san check
        if len(session.cmd.get_object_list()) != 1:
            _LOGGER.warning(
                "more than 1 object is found in current pymol session when doing select. May give absurd results if distance based selection is used. (It will use atoms from other objs)"
            )
        sele_pattern = f"{pymol_obj_name} & ({pattern})"
        sele_name = session.cmd.get_unused_name("enzy_htp_stru_sele")
        session.cmd.select(sele_name, sele_pattern)
        session.cmd.iterate(sele_name, "result.append(ID)", space=locals())  # use ID instead of index here

        return sorted(result)

    def point_mutate(self,
                     pymol_session: pymol2.PyMOL,
                     pos_key: Tuple[str, int],
                     target: str,
                     pymol_obj_name: str,
                     debug: bool = False) -> None:
        """
        Performs a single point mutation on the target WT in the PyMOL session in-place.
        Args:
            pymol_session: the target PyMOL session.
            pos_key: the chain id and residue index.
            target: the target residue name (3 letters).
            pymol_obj_name: the name of the target WT in PyMOL.
            debug: prints all rotamer scores out for debugging purposes.
        Returns:
            None.
        """

        # load pymol wizard mutagenesis function
        pymol_session.cmd.wizard("mutagenesis")
        pymol_session.cmd.do("refresh_wizard")

        # select res idx to mutate to target
        pymol_session.cmd.get_wizard().set_mode(target)
        if pymol_obj_name == "":
            pymol_session.cmd.get_wizard().do_select(f"{pos_key[0]}/{str(pos_key[1])}/")
        else:
            pymol_session.cmd.get_wizard().do_select(f"{pymol_obj_name}//{pos_key[0]}/{str(pos_key[1])}/")

        # prints all rotamers and strains out; also saves each variation to a PDB file in scratch/.
        # can use for debugging purposes
        if debug:
            for i in range(1, pymol_session.cmd.count_states() + 1):
                pymol_session.cmd.wizard("mutagenesis")
                pymol_session.cmd.do("refresh_wizard")
                pymol_session.cmd.get_wizard().set_mode(target)
                if pymol_obj_name == "":
                    pymol_session.cmd.get_wizard().do_select(f"{pos_key[0]}/{str(pos_key[1])}/")
                else:
                    pymol_session.cmd.get_wizard().do_select(f"{pymol_obj_name}//{pos_key[0]}/{str(pos_key[1])}/")

                pymol_session.cmd.get_wizard().do_state(i)
                pymol_session.cmd.frame(i)
                pymol_session.cmd.get_wizard().apply()
                self.export_pdb(pymol_obj_name, pymol_session, tag="rotamer_" + str(i))
            return

        pymol_session.cmd.get_wizard().apply()


    def export_pdb(self, pymol_session: pymol2.PyMOL, pymol_obj_name: str, if_retain_order: bool = True, tag: str = None) -> str:
        """
        Saves a PyMOL object to a PDB file.
        Args:
            pymol_session: the target PyMOL session.
            pymol_obj_name: the name of the target enzyme in PyMOL.
            if_retain_order: if the saving should keep the order of the atoms in the original object.
            tag: the name tag for the saved file.
        Returns:
            The path to the saved PDB file.
        """
        result_dir = eh_config["system.SCRATCH_DIR"]
        fs.safe_mkdir(result_dir)
        if tag is not None:
            pymol_outfile_path = fs.get_valid_temp_name(f"{result_dir}/{pymol_obj_name}_{tag}.pdb")
        else:
            pymol_outfile_path = fs.get_valid_temp_name(f"{result_dir}/{pymol_obj_name}.pdb")

        if if_retain_order:
            pymol_session.cmd.set("retain_order")
        pymol_session.cmd.save(pymol_outfile_path, pymol_obj_name)

        return pymol_outfile_path

    def export_enzy_htp_stru(self, 
                             pymol_session: pymol2.PyMOL,
                             pymol_obj_name: str,        
                             if_retain_order: bool = False) -> Structure:
        """
        Saves a PyMOL object to a Structure object.
        Args:
            pymol_session: the target PyMOL session.
            pymol_obj_name: the name of the target enzyme in PyMOL.
            if_retain_order: if the saving should keep the order of the atoms in the original object.
        Returns:
            A Structure object representing the target enzyme in PyMOL.
        """
        sp = PDBParser()
        pymol_outfile_path = self.export_pdb(pymol_obj_name, pymol_session, if_retain_order=if_retain_order)
        res = sp.get_structure(pymol_outfile_path, allow_multichain_in_atom=True)
        fs.clean_temp_file_n_dir([pymol_outfile_path, eh_config["system.SCRATCH_DIR"]])
        return res
    # == inter-session modular functions == (do not requires a session, will start and close one)
    # pass

    # TODO: go to a parent class
    @property
    def config(self) -> PyMolConfig:
        """Getter for the PyMolConfig() instance belonging to the class."""
        return self.config_

    def de_protonate(self, session: pymol2.PyMOL, molfile: str, new_file: str = None) -> str:
        """Method that deprotonates the structure supplied by molife, saving it to a new file. 
        If no new_file is supplied, original file <name>.<ext> will be saved to <name>_deprotonated.<ext>

        Args:
            session : A pymol2.PyMOL() session to use. 
            molfile: The name of the input file.
            new_file: Where to save the deprotonated structure. Optional.
    
        Returns:
            The name of the deprotonated structure.
        """

        _LOGGER.warning("This function should be used for non-amino acid structures only!")

        self.check_pymol2_installed()

        fs.check_file_exists(molfile)

        if not new_file:
            fpath = Path(molfile)
            new_file = fpath.parent / f"{fpath.stem}_deprotonated{fpath.suffix}"

        session.cmd.delete('all')
        session.cmd.load(molfile)
        session.cmd.remove('hydrogens')
        session.cmd.save(new_file)

        session.cmd.delete('all')

        return new_file

    def rename_atoms(self, template: str, in_file: str) -> str:
        """ """
        _eh_local: Dict[str, Any] = {'template': [], 'in_stru': []}

        self.cmd.delete('all')
        self.cmd.load(template)
        self.cmd.iterate('all', 'template.append((name, elem, ID))', space=_eh_local)
        self.cmd.delete('all')

        self.cmd.load(in_file)
        self.cmd.iterate('all', 'in_stru.append((name, elem, ID))', space=_eh_local)

        if len(_eh_local['template']) != len(_eh_local['in_stru']):
            _LOGGER.error(
                f"Atom number mismatch in template vs in_file structures. ({len(_eh_local['template'])} vs {len(_eh_local['in_stru'])}). Exiting..."
            )
            exit(1)

        for tidx, (tp1, tp2) in enumerate(
                zip(sorted(_eh_local['template'], key=lambda tp: tp[-1]), sorted(_eh_local['in_stru'], key=lambda tp: tp[-1]))):
            if tp1[1] != tp2[1]:
                _LOGGER.error(
                    f"Atom type mismatch at atom number {tidx+1} template vs in_file structures. template: {tp1[1]}  in_stru: {tp2[1]}. Exiting..."
                )
                exit(1)

        for (name, _, ID) in _eh_local['template']:
            self.cmd.alter(f"ID {ID}", f'name="{name}"')

        self.cmd.save(in_file)
        self.cmd.delete('all')

        return in_file

    def collect(self, session: pymol2.PyMOL, molfile: str, variables: List[str], sele: str = 'all', state: int = -1) -> pd.DataFrame:
        """Function that serves as a wrapper to the pymol cmd API's iterate_state() method. Checks if the supplied file exists prior to collecting data.
        Data is returned in a pandas DataFrame() where each row is an atom and each column is a variable.

        Args:
            session : A pymol2.PyMOL() session to use. 
            molfile: Name of the structure to sample. 
            variables: A list() of str()'s with variable names to pull from in the structure.
            sele: The pymol selection string to use. Defaults to 'all'
            state: The state to use. Defaults to -1 (current state).

        Returns:
            A pandas DataFrame() where each row is an atom and each column is a variable.

        """
        self.check_pymol2_installed()

        check_var_type(variables, list)
        for vv in variables:
            check_var_type(vv, str)

        if molfile != 'memory':
            fs.check_file_exists(molfile)
            session.cmd.delete('all')
            session.cmd.load(molfile)

        _eh_local: Dict[str, Any] = {'data': []}

        iter_stmt: str = f"data.append(({', '.join(variables)}))"
        session.cmd.iterate_state(state, sele, iter_stmt, space=_eh_local)

        result = pd.DataFrame(data=_eh_local['data'], columns=variables)

        return result

    def general_cmd(
        self,
        session: pymol2.PyMOL,
        args: List[Tuple],
    ) -> List[Any]:
        """Executes a series of commands through the PyMOL/PyMOL2 python module in use. Takes input as a list of Tuple()'s
        where the first item in each tuple is a string specifying the function to use and the rest of the items are the
        arguments for that function.
        TODO(CJ): add examples
    
        Args:
            session : A pymol2.PyMOL() session to use.
            args: A list() of tuple()'s of items where the first string names the pymol API command to use and the rest are arguments for that function.

        Returns:
            Each result for the respective command output.
        """

        self.check_pymol2_installed()

        result: List[Any] = list()

        for cmd_set in args:
            cmd_name = cmd_set[0]

            #if cmd_name not in self.available_cmds_:
               # _LOGGER.error(f"The command '{cmd_name}' is not supported in this version of pymol. Exiting...")
                #exit(1)
            
            cmd_args = list()

            if len(cmd_set) > 1:
                cmd_args = list(cmd_set[1:])

            cmd_str: str = f"{cmd_name}({','.join(map(str, cmd_args))})"
            try:
                fxn = getattr(session.cmd, cmd_name)
                if len(cmd_set) > 1:
                    result.append(fxn(*cmd_args))
                else:
                    result.append(fxn())
            except Exception as e:
                _LOGGER.error(f"{e}")
                _LOGGER.error(f"PyMOL function call '{cmd_str}' resuled in an error. Exiting...")
                exit(1)

        return result

    def get_charge(self, session: pymol2.PyMOL, fname: str, sele: str = '(all)') -> int:
        """Method that gets the formal charge for the specified sele in the specified file. File must be a supported file type as listed in
        PyMOLConfig. Checks if file exists and is supported type. If either
        are not true then an error is logged and the program exits. DOES
        NOT check if the sele is valid.
        Args:
            session : A pymol2.PyMOL() session to use. 
            fname: The str() path to the file in question.
            sele: The str() specifying the atom selection in PyMOL synatx. Default is '(all)'.

        Returns:
            The formal charge of the sele as an int().
        """

        self.check_pymol2_installed()

        if not Path(fname).exists():
            _LOGGER.error(f"The supplied file '{fname}' does not exist. Exiting...")
            exit(1)

        if not self.supported_file_type(fname):
            _LOGGER.error(f"The supplied file '{fname}' has an unsupported extension. Exiting...")
            exit(1)

        _eh_local: Dict[str, Any] = {'fc': []}

        session.cmd.delete('all')
        session.cmd.load(fname)
        session.cmd.iterate(sele, 'fc.append(formal_charge)', space=_eh_local)
        session.cmd.delete('all')

        return sum(_eh_local['fc'])

    def get_sequence(self, session: pymol2.PyMOL, fname: str, sele: str = '(all)') -> str:
        """Method that gets the sequence for the specified sele in the
        specified file. File must be a supported file type as listed in
        PyMOLConfig. Checks if file exists and is supported type. If either
        are not true then an error is logged and the program exits. DOES
        NOT check if the sele is valid. Note that non-canonical residues will
        be represented by a '?' in the sequence.
        Args:
            session : A pymol2.PyMOL() session to use. 
            fname: The str() path to the file in question.
            sele: The str() specifying the atom selection in PyMOL synatx. Default is '(all)'.

        Returns:
            A str() with the sequence of the sele as it would appear in a .fasta file..
        """

        self.check_pymol2_installed()

        if not Path(fname).exists():
            _LOGGER.error(f"The supplied file '{fname}' does not exist. Exiting...")
            exit(1)

        if not self.supported_file_type(fname):
            _LOGGER.error(f"The supplied file '{fname}' has an unsupported extension. Exiting...")
            exit(1)

        session.cmd.delete('all')
        session.cmd.load(fname)
        lines: List[str] = session.cmd.get_fastastr(sele).splitlines()
        session.cmd.delete('all')

        lines = list(filter(lambda ll: ll[0] != '>', lines))
        return ''.join(lines)


    def create_cluster(self,session, fname:str, sele_str:str,outfile:str=None,cap_strategy:str='H', work_dir:str=None) -> str:
        """ 
        """
    
        if work_dir is None:
            work_dir = config['system.WORK_DIR']
    
        fs.check_file_exists(fname)
    
        if outfile is None:
            temp_path = Path(fname)
            outfile:str = f"{work_dir}/{temp_path.stem}_cluster{temp_path.suffix}"
    
        #TODO(CJ): add file check for fname
        obj_name:str='__eh_cluster'
        self.general_cmd(session, [('load', fname), ('select', sele_str),('create',obj_name,sele_str)])
    
        df:pd.DataFrame=self.collect(session, 'memory', "chain resi resn name".split())
    
        args = list()
        for i,row in df.iterrows():
            if row['name'] not in "N C".split():
                continue
            
            if not row.resn in THREE_LETTER_AA_MAPPER:
                continue
            
            args.append(
                ('h_add', f"chain {row.chain} and resi {row.resi} and resn {row.resn} and name {row['name']}")
            )
    
        args.append(("save", outfile, obj_name))
    
        self.general_cmd(session, args )
    
        if cap_strategy == 'H':
            return outfile
    
        if cap_strategy == 'CH3':
            
            args = []
            orig = set(
                zip(df.chain,df.resi,df.resn,df['name'])
            )
            updated:pd.DataFrame=self.collect(session, 'memory', "chain resi resn name".split())
            new = set(
                zip(updated.chain,updated.resi,updated.resn,updated['name'])
            )
    
            for tp in new:
                if tp in orig:
                    continue
                if tp[3] == 'H01':
                    new_name='C21'
                elif tp[3] == 'H02':
                    new_name='C22'
                else:
                    assert False
                args.extend([
                    ('alter', f"chain {tp[0]} and resi {tp[1]} and resn {tp[2]} and name {tp[3]}", "elem='C'"),
                    ('alter', f"chain {tp[0]} and resi {tp[1]} and resn {tp[2]} and name {tp[3]}", f"name='{new_name}'"),
                ])
    
            args.extend([
                ('valence', 'guess', 'name C21 or name C22'),
                ('h_add','name C21'),
                ('h_add','name C22'),
                ("save", outfile, obj_name)
            ])
            self.general_cmd(session, args )
    
        return outfile    


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
