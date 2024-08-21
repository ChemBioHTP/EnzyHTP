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
import numpy as np
import pandas as pd

#TODO(CJ): add something to remove "PyMOL not running. Entering library mode (experimental" message on pymol running
from enzy_htp import config as eh_config
from enzy_htp.core import env_manager as em
from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER, check_var_type
from enzy_htp._config.pymol_config import PyMolConfig, default_pymol_config
from enzy_htp.structure import Structure, PDBParser, Ligand
import enzy_htp.chemical as chem
from enzy_htp.structure.structure_ensemble import StructureEnsemble

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
            _LOGGER.error(f"pymol2 is NOT installed. Use 'conda install -c conda-forge -y -q pymol-open-source'. Exiting...")
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

        self.general_cmd(session, [('delete', 'all')])

        try:
            self.general_cmd(session, [('load', file_1)])

            if split_states:
                self.general_cmd(session, [('split_states', 'all')])
                result = Path(result)
                parent_dir = result.parent
                result_stem = result.stem
                result = []

                for oidx, oname in enumerate(self.general_cmd(session, [('get_object_list')])[0]):

                    if self.general_cmd(session, [('count_states', oname)])[0] != 1:
                        session.cmd.delete(oname)
                        continue

                    fname = str(parent_dir / f"{result_stem}_{oidx}{new_ext}")
                    self.general_cmd(session, [('save', fname, oname)])
                    result.append(fname)
            else:
                self.general_cmd(session, [('save', result), ('delete', 'all')])

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
        """convert enzy_htp.Structure into a pymol object in a pymol2.PyMOL() session. Performs 
        operation without writing to disk.
        
        Returns:
            (pymol_obj_name, session) since the name is only valid in the session
            
        Note: pymol wont reset residue idx or chain names but will reset atom index from 1"""

        self.check_pymol2_installed()
        
        pdb_str:str = PDBParser().get_file_str(stru, if_renumber=False)
        
        pymol_obj_name:str = self.general_cmd(session, [('get_unused_name',)])[0] 
        
        self.general_cmd(session, [('read_pdbstr', pdb_str, pymol_obj_name)])


        return (pymol_obj_name, session)
    
    def load_enzy_htp_stru_esm(self, session: pymol2.PyMOL, stru_esm: StructureEnsemble) -> Tuple[str, pymol2.PyMOL]:
        """Convert enzy_htp.structure.StructureEnsemble object into a pymol object in a pymol2.PyMOL() session.
        
        Args:
            session (pymol2.PyMOL): Current pymol2.PyMOL() session instance.
            stru_esm (StructureEnsemble): A StructureEnsemble instance.

        Returns:
            pymol_obj_name (str): The name of the loaded PyMOL object.
            session (pymol2.PyMOL): Current pymol2.PyMOL() session instance.
        """
        self.check_pymol2_installed()
        prmtop_str = stru_esm._topology

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

    def export_enzy_htp_stru(self, pymol_obj_name: str,
                            pymol_session: pymol2.PyMOL,
                            if_retain_order: bool = False,
                            if_fix_naming: bool = False) -> Structure:
        """
        Saves a PyMOL object to a Structure object.
        Args:
            pymol_session: the target PyMOL session.
            pymol_obj_name: the name of the target enzyme in PyMOL.
            if_retain_order: if the saving should keep the order of the atoms in the original object.
                             Fixes PyMOL's scrambling of the atom order
            if_fix_naming: if the PyMOL naming of atoms should be updated to match the Structure convention.
                           Allows for atom naming consistency (https://github.com/ChemBioHTP/EnzyHTP/issues/117)
        Returns:
            A Structure object representing the target enzyme in PyMOL.
        """
        sp = PDBParser()
        pymol_outfile_path = self.export_pdb(pymol_session, pymol_obj_name,
                                             if_retain_order=if_retain_order)
        res = sp.get_structure(pymol_outfile_path, allow_multichain_in_atom=True)

        if if_fix_naming:
            self.fix_pymol_naming(res)

        fs.clean_temp_file_n_dir([pymol_outfile_path, eh_config["system.SCRATCH_DIR"]])
        return res


    def fix_pymol_naming(self, stru: Structure) -> None:
        """Fixes PyMOL's naming of atoms in-place using the PYMOL_TO_ATOM_MAPPER in pymol_config.py.
        Args:
            stru: the Structure that contains the wrong atom names (if any).
        Returns:
            Nothing.
        """
        for chain in stru.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    correct_name = PyMolConfig().get_canonical_atom_name(residue.name, atom.name)
                    if correct_name:
                        atom.name = correct_name

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
                if cmd_name == 'alter' and len(cmd_set) == 3:
                    #TODO(CJ): add some more checking here
                    result.append(session.cmd.alter(
                        cmd_set[1], cmd_set[2]
                    ))
                elif len(cmd_set) > 1:
                    
                    result.append(fxn(*cmd_args))
                else:
                    result.append(fxn())
            except Exception as e:
                _LOGGER.error(f"{e}")
                _LOGGER.error(f"PyMOL function call '{cmd_str}' resulted in an error. Exiting...")
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

    def remove_ligand_bonding(self, session: pymol2.PyMOL, clash_cutoff:float=2.0) -> None: #TODO(CJ): bad name here
        """Given a session with some already present molecule, remove the bonding in between residues presumably due
        to clashes. Essentially loops through all atoms and unbonds if they are within the specified clash_cutoff.

        Args:
            session: The pymol session all molecular information is stored inside of.
            clash_cutoff: How close can two atoms be before they are bonded? Default is 2.0 angstroms.

        Returns:
            Nothing.
        """

        self.general_cmd(session, [('set', 'retain_order', 1)])       #NOTE(CJ): this solves so many problems
        #TODO(CJ): need to add 
        def dist( p1, p2 ):
            return np.sqrt(np.sum((p1-p2)**2))
        #TODO(CJ): documentation
        df: pd.DataFrame = self.collect(session, 'memory', "chain resi resn name x y z".split(), "not polymer.protein")
        df['point'] = df.apply(lambda row: np.array([row.x, row.y, row.z]) ,axis = 1)
        non_aa = set(list(zip(df.chain, df.resi, df.resn)))
        sele_str = " or ".join(map(lambda pr: f"(chain {pr[0]} and resi {pr[1]} and resn {pr[2]})", non_aa))

        df2 = self.collect(session, 'memory', 'chain resi resn name x y z'.split(), f'polymer.protein within 3 of {sele_str}')
        df2['point'] = df2.apply(lambda row: np.array([row.x, row.y, row.z]) ,axis = 1)
       
        args = list()
        for i, row in df.iterrows():
            for i2, row2 in df2.iterrows():
                if dist(row.point, row2.point) <= clash_cutoff:
                    args.append(('unbond',
                        f"chain {row.chain} and resi {row.resi} and resn {row.resn}",
                        f"chain {row2.chain} and resi {row2.resi} and resn {row2.resn}",
                    ))
        self.general_cmd(session, args)


    def create_cluster(self, session, fname: str, sele_str: str, outfile: str = None, cap_strategy: str = 'H', work_dir: str = None) -> str:
        """TODO(CJ)"""

        if work_dir is None:
            work_dir = eh_config['system.WORK_DIR']

        fs.check_file_exists(fname)

        if outfile is None:
            temp_path = Path(fname)
            outfile: str = f"{work_dir}/{temp_path.stem}_cluster{temp_path.suffix}"

        #TODO(CJ): add file check for fname
        obj_name: str = '__eh_cluster'
        self.general_cmd(session, [('delete', 'all'), ('load', fname), ('select', sele_str), ('create', obj_name, sele_str),
                                   ('delete', Path(fname).stem)])

        self._remove_ligand_bonding(session)

        df: pd.DataFrame = self.collect(session, 'memory', "chain resi resn name".split())

        args = list()
        for i, row in df.iterrows():
            if row['name'] not in "N C".split():
                continue

            if not row.resn in chem.THREE_LETTER_AA_MAPPER:
                continue

            args.extend([('valence', 'guess', f"chain {row.chain} and resi {row.resi} and resn {row.resn} and name {row['name']}"),
                         ('h_add', f"chain {row.chain} and resi {row.resi} and resn {row.resn} and name {row['name']}")])

        args.append(("save", outfile, obj_name))

        self.general_cmd(session, args)

        if cap_strategy == 'H':
            return outfile

        if cap_strategy == 'CH3':

            args = []
            orig = set(zip(df.chain, df.resi, df.resn, df['name']))
            updated: pd.DataFrame = self.collect(session, 'memory', "chain resi resn name".split())
            new = set(zip(updated.chain, updated.resi, updated.resn, updated['name']))

            for (cname, res_num, res_name, aname) in filter(lambda x: x not in orig, new):
                if aname == 'H01':
                    new_name = 'C21'
                elif aname == 'H02':
                    new_name = 'C22'
                else:
                    #TODO(CJ): better error message
                    self.general_cmd(session, [('save', '_____mess_up.pdb')])
                    assert False, (cname, res_num, res_name, aname)
                args.extend([
                    ('alter', f"chain {cname} and resi {res_num} and resn {res_name} and name {aname}", "elem='C'"),
                    ('alter', f"chain {cname} and resi {res_num} and resn {res_name} and name {aname}", f"name='{new_name}'"),
                ])

            self.general_cmd(session, args)
            args = [('valence', 'guess', 'name C21 or name C22'), ('h_add', 'name C21'), ('h_add', 'name C22'), ("save", outfile, obj_name)]
            self.general_cmd(session, args)

        return outfile


    def center_of_mass(self, session, sele:str='all', no_hydrogens:bool=True):
        """TODO(CJ)
        Args:
            session:
            sele:
            no_hydrogens:

        Returns:
            The specified center of mass.
        """
        if no_hydrogens:
            sele += " and (not elem H)"
            
        df = self.collect(session, 'memory', "x y z".split(), sele=sele)
        return np.mean(np.array([df.x.to_numpy(), df.y.to_numpy(), df.z.to_numpy()]),axis=1)

    def fetch(self, code: str, out_dir: str = None) -> str:
        """Given a 

        Args:
            code:
            out_dir

        Returns:
            The path to the
        """
        outfile = f"{code.upper()}.cif"
        if len(code) == 3:
            url: str = f"{self.config_.LIGAND_STEM}/{outfile}"
        elif len(code) == 4:
            url: str = f"{self.config_.STRUCTURE_STEM}/{outfile}"
        else:
            assert False

        self.env_manager_.run_command(self.config_.WGET, [url])

        if out_dir is not None:
            outfile = fs.safe_mv(outfile, f"{out_dir}/")
        #TODO(CJ): check if the file is downloaded
        return outfile


    

    def get_residue_list(self, session, stru, sele_str:str='all', work_dir:str=None) -> List[Tuple[str,int]]:
        #TODO(CJ): add this documentation + type hinting
        if work_dir is None:
            work_dir = './'

        temp_file:str = f"{work_dir}/__temp_pymol.pdb"

        _parser = PDBParser()
        _parser.save_structure(temp_file, stru)

        df = self.collect(session, temp_file, "chain resi".split(), sele=sele_str)
        
        fs.safe_rm( temp_file )
        result = list()
        result_set = set()

        for i, row in df.iterrows():
            new = (row['chain'], int(row['resi']))
            if new not in result_set:
                result.append( new )
                result_set.add( new )

        return result 

    
    def get_atom_mask(self, session, stru:Structure, sele_str:str=None, work_dir:str=None) -> List[bool]:
        #TODO(CJ): this stuff

        
        if work_dir is None:
            work_dir = './'

        temp_file:str = f"{work_dir}/__temp_pymol.pdb"
        _parser = PDBParser()
        _parser.save_structure(temp_file, stru)

        session = self.new_session()

        df:pd.DataFrame = self.collect(session, temp_file, "chain resi name".split(), sele=sele_str)
        
        fs.safe_rm( temp_file )

        sele_set = set()

        for i, row in df.iterrows():
            sele_set.add(f"{row['chain']}.{row['resi']}.{row['name']}")

        mask:List[bool] = list()

        for atom in stru.atoms:
            mask.append( atom.key in sele_set )            

        return mask

    def get_rmsd(self, structure_ensemble: StructureEnsemble, mask_pattern: str = str()) -> float:
        """Calculates the Root Mean Square Deviation (RMSD) for a structure ensemble. Ligands, solvents and ions are not included.
        The lower the RMSD value is, the more stable the structure is.

        Args:
            structure_ensemble (StructureEnsemble): A collection of different geometries of the same enzyme structure.
            mask_pattern (str, optional): A pymol-formatted selection string which defines the region for calculating RMSD value.

        Returns:
            The calculated RMSD value as a float.
        """
        # TODO: By residue.
        with OpenPyMolSession(self) as pms:
            ref_stru_name, _ = self.load_enzy_htp_stru(pms, structure_ensemble.structure_0)
            rmsd_list = list()
            for frame_stru in structure_ensemble.structures(remove_solvent=True):
                frame_stru_name, _ = self.load_enzy_htp_stru(pms, stru=frame_stru)
                results: List[Any] = self.general_cmd(pms, [
                    ('remove', 'solvent'),
                    ('remove', 'inorganic'),
                    ('remove', 'organic'),
                    ('sele', mask_pattern),

                    ('align', frame_stru_name, ref_stru_name),
                    # This returns a tuple with 7 items:
                    # 1. RMSD after refinement
                    # 2. Number of aligned atoms after refinement
                    # 3. Number of refinement cycles
                    # 4. RMSD before refinement
                    # 5. Number of aligned atoms before refinement
                    # 6. Raw alignment score
                    # 7. Number of residues aligned
                    # https://pymolwiki.org/index.php/Align

                    ('remove', frame_stru_name),
                ])
                rmsd_list.append(results[-1][0])
                continue
            return sum(rmsd_list) / len(rmsd_list)
        
    def get_spi(self, stru: Structure, ligand:Ligand, pocket_sele:str) -> float:
        """Calculates substrate positioning index (SPI) for a given Ligand in a given Structure. SPI roughly corresponds to the
        ratio of the ligand's solvent accessible surface area (SASA) dived by protein binding pocket SASA. The citation for 
        this paper is here: DOI:https://doi.org/10.1021/acs.jpclett.3c02444.

        Args:
            stru: The Structure containing the Ligand and active site.
            ligand: The actual Ligand we are calculating SPI for. 
            pocket_sele: A str describing the active site of the protein in pymol format.
        
        Returns:
            The calculated SPI as a float.

        """
        (lig_chain, lig_idx) = ligand.key()
        with OpenPyMolSession(self) as pms:
            self.load_enzy_htp_stru(pms, stru )
            results:List[Any] = self.general_cmd(pms,[
                ('set', 'dot_solvent', 1),
                ('create', 'ligand', f'chain {lig_chain} and resi {lig_idx} and not solvent'),
                ('create', 'protein', f'not (chain {lig_chain} and resi {lig_idx}) and not solvent'),
                ('get_area', 'ligand'),
                ('get_area', f'protein and ({pocket_sele})')
            ])
            return results[-2] / results[-1]

class OpenPyMolSession:
    """a context manager that open a pymol session once enter and close once exit"""

    def __init__(self, pymol_interface: PyMolInterface) -> None:
        self.interface = pymol_interface

    def __enter__(self):
        """open a pymol session once enter"""
        self.session = self.interface.new_session()
        return self.session

    def __exit__(self, exc_type, exc_val, exc_tb):
        """open a pymol session once enter"""
        self.session.stop()
