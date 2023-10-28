"""Defines an AmberInterface class that serves as a bridge for enzy_htp to utilize AmberMD software. Uses the AmberConfig class
found in enzy_htp/molecular_mechanics/amber_config.py. Supported operations include mutation with tLEaP, MolDynParameterizer for MD,
and MolDynStep for modular MD steps that can be minimization, heating, constant pressure production, or constant pressure
equilibration.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-02
"""
import glob
import os
import re
import shutil
from pathlib import Path
from subprocess import CalledProcessError
from typing import List, Tuple, Union, Dict, Any

from .base_interface import BaseInterface
from .handle_types import MolDynParameterizer, MolDynParameter, MolDynStep
from .ncaa_library import search_ncaa_parm_file
from .gaussian_interface import gaussian_interface

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.core import env_manager as em
from enzy_htp.core.exception import UnsupportedMethod, tLEaPError
from enzy_htp._config.amber_config import AmberConfig, default_amber_config
from enzy_htp.structure.structure_io import pdb_io
from enzy_htp.structure import (
    Structure,
    Ligand,
    MetalUnit,
    ModifiedResidue,
    NonCanonicalBase,)
from enzy_htp import config as eh_config

class AmberParameter(MolDynParameter):
    """the Amber format MD parameter. Enforce the Amber parameter format
    in AmberMDStep.
    Attribute:
        inpcrd: the path of the .inpcrd input coordinate file
        prmtop: the path of the .prmtop parameter topology file"""

    def __init__(self, inpcrd_path: str, prmtop_path: str):
        self._inpcrd = inpcrd_path
        self._prmtop = prmtop_path

    #region == property ==
    @property
    def engine(self) -> str:
        return "Amber"

    def get_solvated_structure(self) -> Structure:
        """get the solvated structure corresponding to the parameters"""
        pass # TODO do we really need this?

    @property
    def file_list(self) -> List[str]:
        """return a list of files that composes the parameter"""
        return [self._inpcrd, self._prmtop]
    #endregion

    #region == checker ==
    def is_valid(self) -> bool:
        """check whether the parameter files represented by {self} is valid"""
        result = 1
        # file exist
        result *= Path(self._inpcrd).exists()
        result *= Path(self._prmtop).exists()
        # file size not zero
        result *= os.path.getsize(self._inpcrd) != 0
        result *= os.path.getsize(self._prmtop) != 0
        # TODO add upon need
        return bool(result)

    #endregion


class AmberParameterizer(MolDynParameterizer):
    """the MD parameterizer for Amber.
    Constructer:
        AmberInterface.build_md_parameterizer()
    Attributes: (configuration of the parameterization process)
        charge_method
        resp_engine
        resp_lvl_of_theory
        ncaa_param_lib_path
        force_renew_ncaa_parameter
        ncaa_net_charge_engine
        ncaa_net_charge_ph
        solvate_box_type
        solvate_box_size
        gb_radii
        parameterizer_temp_dir
        additional_tleap_lines"""

    def __init__(
            self,
            force_fields: List[str],
            charge_method: str,
            resp_engine: str,
            resp_lvl_of_theory: str,
            ncaa_param_lib_path: str,
            force_renew_ncaa_parameter: bool,
            ncaa_net_charge_engine: str,
            ncaa_net_charge_ph: float,
            solvate_box_type: str,
            solvate_box_size: float,
            gb_radii: str,
            parameterizer_temp_dir: str,
            additional_tleap_lines: List[str],) -> None:
        self.force_fields = force_fields
        self.charge_method = charge_method
        self.resp_engine = resp_engine
        self.resp_lvl_of_theory = resp_lvl_of_theory
        self.ncaa_param_lib_path = ncaa_param_lib_path
        self.force_renew_ncaa_parameter = force_renew_ncaa_parameter
        self.ncaa_net_charge_engine = ncaa_net_charge_engine
        self.ncaa_net_charge_ph = ncaa_net_charge_ph
        self.solvate_box_type = solvate_box_type
        self.solvate_box_size = solvate_box_size
        self.gb_radii = gb_radii
        self.parameterizer_temp_dir = parameterizer_temp_dir
        self.additional_tleap_lines = additional_tleap_lines

    @property
    def engine(self) -> str:
        return "Amber"

    def run(self, stru: Structure) -> AmberParameter:
        """the parameterizer convert stru to amber parameter (inpcrd, prmtop)"""
        # 0. set up paths
        result_inpcrd = fs.get_valid_temp_name(f"{self.parameterizer_temp_dir}/amber_parm.inpcrd")
        result_prmtop = fs.get_valid_temp_name(f"{self.parameterizer_temp_dir}/amber_parm.prmtop")
        fs.safe_mkdir(self.parameterizer_temp_dir)

        # 1. check stru diversity
        diversity = stru.chemical_diversity
        _LOGGER.debug(f"diversity: {diversity}")

        # 2. extract and parameterize each special component
        ligand_parms = {}
        maa_parms = {}
        metalcenter_parms = {}
        gaff_type = self._check_gaff_type()
        if "ligand" in diversity:
            for lig in stru.ligands:
                if lig.name not in ligand_parms: # avoid repeating calculation
                    ligand_parms[lig.name] = self._parameterize_ligand(lig, gaff_type)
        if "modified_residue" in diversity:
            for maa in stru.modified_residue:
                if maa.name not in maa_parms:
                    maa_parms[maa.name] = self._parameterize_modified_res(maa, gaff_type)
        if "metalcenter" in diversity:
            _LOGGER.warning(
                "Support for paramization of metalcenter is not ready yet."
                " You can use AmberParameterizer().additional_tleap_lines"
                " & change relevant residue names to apply pre-made parameters"
                " for metal centers.")
            # all_used_residue_names = []
            # all_used_atom_types = []
            # for metal in stru.metalcenters:
            #     parms, used_residue_names, used_atom_types = self._parameterize_metalcenter(
            #         metal, ligand_parms, maa_parms)
            #     metalcenter_parms.update(parms)
            #     all_used_residue_names.extend(used_residue_names)
            #     all_used_atom_types.extend(used_atom_types)

        # 3. write the combining tleap.in
        tleap_content, temp_dry_pdb = self._write_combining_tleap_input(
                                            stru,
                                            ligand_parms,
                                            maa_parms,
                                            metalcenter_parms,
                                            result_inpcrd,
                                            result_prmtop,)

        # 4. run tleap
        amber_interface.run_tleap(tleap_content)

        # 5. clean up
        fs.clean_temp_file_n_dir([
            temp_dry_pdb,
            self.parameterizer_temp_dir,
        ])

        return AmberParameter(result_inpcrd, result_prmtop)

    def _parameterize_ligand(self, lig: Ligand, gaff_type: str) -> Tuple[str, List[str]]:
        """parameterize ligand for AmberMD, use ncaa_param_lib_path for customized
        parameters. Multiplicity and charge information can be set in Ligand objects.
        TODO prob add a Structure method for batch assigning it
        Returns:
            mol_desc_path, [frcmod_path, ...]"""
        # san check
        if not gaff_type:
            _LOGGER.error("The structure contains non-canonical residue"
                          " (lig or maa) but GAFF/GAFF2 is not used!"
                          f" Check you force_fields. (current: {self.force_fields})")
            raise ValueError
        # init
        fs.safe_mkdir(self.ncaa_param_lib_path)
        target_method = f"{self.charge_method}-{gaff_type}"

        # 0. search parm lib
        mol_desc_path, frcmod_path_list = search_ncaa_parm_file(lig,
                                            target_method=target_method,
                                            ncaa_lib_path=self.ncaa_param_lib_path)

        if mol_desc_path:
            if frcmod_path_list:
                return mol_desc_path, frcmod_path_list
        else:
            # 1. generate mol2 if not found
            mol_desc_path = f"{self.ncaa_param_lib_path}/{lig.name}_{target_method}.mol2" # the search ensured no existing file named this
            amber_interface.antechamber_ncaa_to_moldesc(ncaa=lig,
                                                        out_path=mol_desc_path,
                                                        gaff_type=gaff_type)

        # 2. run parmchk2 on the PDB
        # (this also runs when mol_desc exsit but not frcmod)
        frcmod_path = f"{self.ncaa_param_lib_path}/{lig.name}_{target_method}.frcmod"
        amber_interface.run_parmchk2(in_file=mol_desc_path,
                                     out_file=frcmod_path,
                                     gaff_type=gaff_type)

        return mol_desc_path, [frcmod_path]

    def _parameterize_modified_res(self, maa: ModifiedResidue, gaff_type: str) -> Tuple[str, List[str]]:
        """parameterize modified residues for AmberMD, use ncaa_param_lib_path for customized
        parameters. Multiplicity and charge information can be set in ModifiedResidue objects."""
        fs.safe_mkdir(self.ncaa_param_lib_path)
        # 0. search parm lib

        # 1. make maa PDB

        # 2. run antechamber on the PDB get ac
        ac_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.ac")
        # 2.1 fix the wrong atom type given by antechamber

        # 3. run prepgen on ac & mc get prepin
        prepin_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.prepin")

        # 4. run antechamber on prepin get mol2
        mol2_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.mol2")

        # 5. run parmchk2 twice on prepin get frcmod & frcmod2
        frcmod_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.frcmod")
        frcmod2_path = fs.get_valid_temp_name(
            f"{self.ncaa_param_lib_path}/{maa.name}.frcmod2")

        return mol2_path, [frcmod_path, frcmod2_path] # TODO make sure whether mol2 works or do we even need it?

    def _parameterize_metalcenter(self, metal: MetalUnit,
                                  ligand_parms: Dict[str, Tuple[str, List[str]]],
                                  maa_parms: Dict[str, Tuple[str, List[str]]],) -> List :
        """TODO maybe rewrote MCPB.py"""
        mcpb_path = f"{self.ncaa_param_lib_path}/mcpb"
        fs.safe_mkdir(mcpb_path)
        # 1. run MCPB.py 1st step - get involved residues from "small", new res name
        mcpbin_path = fs.get_valid_temp_name(
            f"{mcpb_path}/{metal.name}.mcpbin")
        # 2. search parm lib TODO: prob by coordinate

        # 3. run gaussian calculation

        # 4. run MCPB 2/3/4 step get mol2 files, frcmod file, new pdb, new atom types, bond lines
        new_residue_names = []
        new_atom_types = []

        # 5. store mol2, frcmod, bond lines in parm_dict
        parm_dict = {}

        return parm_dict, new_residue_names, new_atom_types

    def _write_combining_tleap_input(self,
            stru: Structure,
            ligand_parms: Dict,
            maa_parms: Dict,
            metalcenter_parms: Dict,
            result_inpcrd: str,
            result_prmtop: str,) -> Tuple[str, str]:
        """combine mol_desc and parm file of each noncanonical parts and make the content
        of the input file for tleap.
        Returns:
            (tleap.in content, temp_pdb_path)"""
        # init file path
        temp_pdb_path = fs.get_valid_temp_name(f"{self.parameterizer_temp_dir}/tleap_combine.pdb")
        # make content
        lines = []
        for ff in self.force_fields:
            lines.append(f"source {ff}")

        # NCAA parts

        # # add atom types TODO
        # for xxx, new_atom_types in metalcenter_parms.values():
        #     pass

        # ligand
        for ncaa_name, (mol_desc_path, frcmod_path_list) in ligand_parms.items():
            lines.extend(
                self._make_ncaa_tleap_lines(ncaa_name, mol_desc_path, frcmod_path_list))
        # MAA
        for ncaa_name, (mol_desc_path, frcmod_path_list) in maa_parms.items():
            lines.extend(
                self._make_ncaa_tleap_lines(ncaa_name, mol_desc_path, frcmod_path_list))
            
        # load PDB
        pdb_io.PDBParser().save_structure(outfile=temp_pdb_path,
                                          stru=stru,
                                          if_renumber=False)
        lines.extend([
            f"a = loadpdb {temp_pdb_path}",
            "center a",
            "addions a Na+ 0",
            "addions a Cl- 0",
            f"solvate{self.solvate_box_type} a TIP3PBOX {self.solvate_box_size}",
            f"saveamberparm a {result_prmtop} {result_inpcrd}",
            "quit",
        ])

        result = "\n".join(lines)

        return result, temp_pdb_path

    def _make_ncaa_tleap_lines(self,
            ncaa_name: str,
            mol_desc_path: str,
            frcmod_path_list: List[str]) -> List[str]:
        """make ncaa moldesc and params lines in tleap.in"""
        result = []
        # params
        for frcmod in frcmod_path_list:
            result.append(f"loadAmberParams {frcmod}")
        # mol desc
        if fs.get_file_ext(mol_desc_path) in [".prepin", ".prepi"]:
            result.append(f"loadAmberPrep {mol_desc_path}")
        elif fs.get_file_ext(mol_desc_path) in [".mol2"]:
            result.append(f"{ncaa_name} = loadmol2 {mol_desc_path}")
        else:
            _LOGGER.error(f"Got file with wrong extension {mol_desc_path}")
            raise ValueError
        return result

    def _check_gaff_type(self) -> str:
        """check the GAFF type used for parameterization.
        Return None if non of them are used."""
        for ff in self.force_fields:
            pattern_match = re.search(r"[Gg][Aa][Ff][Ff]2?", ff)
            if pattern_match:
                return pattern_match.group().upper()
        return None


class AmberMDStep(MolDynStep):
    """the modular MD step of Amber.
    Attributes: (necessary information of the modular MD step)
        pass"""

    def __init__(self) -> None:
        pass

    @property
    def engine(self) -> str:
        return "Amber"


class AmberInterface(BaseInterface):
    """Class that provides a direct inteface for enzy_htp to utilize AmberMD software.
    Main supported operations:
    MolDynStep: a modular MD step.
    MolDynParameterizer: parameterizer for MD.
    run_tleap: tleap interface and a series of tleap based operations like tleap_clean_up_stru.
    TODO: add more
    * Users should use this class as the only way to interact with any functionality
      in Amber or associated tools like tleap or cpptraj.

    Attributes:
        config_	: The AmberConfig() class which provides settings for both running Amber and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() that interface with the shell and ensure all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is compatible with the object itself.
    """

    AMBER_FILE_FORMAT_MAPPER = {
        ".prepin" : "prepi",
    }
    """the mapper for {ext:format} for file formats defined/used by amber
     that are different than the extension TODO add upon need"""

    SUPPORTED_CHARGE_METHOD_MAPPER = {
        "AM1BCC" : "bcc", "bcc" : "bcc",
        "RESP" : "resp", "resp" : "resp",
        "rc" : "rc",
    }
    """dictionary that maps keywords to charge method that current supported in run_antechamber()"""

    def __init__(self, parent, config: AmberConfig = None) -> None:
        """Simplistic constructor that optionally takes an AmberConfig object as its only argument.
        Calls parent class."""
        super().__init__(parent, config, default_amber_config)

    # == interface general ==
    def display_config(self) -> None:
        """Prints all settings for the object's AmberConfig() inteface to stdout using AmberConfig.display()."""
        self.config_.display()

    # == general amber app interface ==
    def get_file_format(self, fname: str) -> str:
        """determine the file type give file path"""
        ext = fs.get_file_ext(fname)
        if "frcmod" in ext:
            return "frcmod"
        return self.AMBER_FILE_FORMAT_MAPPER.get(ext, ext[1:])

    # -- tleap --
    def run_tleap(
        self,
        tleap_in_str: str,
        # cmd based
        if_ignore_start_up: bool = True,
        additional_search_path: List[str] = None,
        tleap_out_path: str = None,
    ) -> None:
        """the python wrapper of running tleap
        Args:
            tleap_in_str:
                the str content of leap.in file
            if_ignore_start_up:
                if adding the '-s' flag to ignore leaprc startup file.
            additional_search_path:
                list of addition search path of leaprc files. Each path
                in the list will be applied using a '-I' flag.
            tleap_out_path:
                file path for stdout of the tleap command. the _LOGGER level
                determines if delete the file.

        NOTE: run_tleap API should not handle the index alignment since it do
        not carry information of the input pdb"""
        temp_path_list = []
        # init file paths (tleap_in_path, tleap_out_path)
        fs.safe_mkdir(eh_config["system.SCRATCH_DIR"])
        tleap_in_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/tleap.in")
        temp_path_list.extend([eh_config["system.SCRATCH_DIR"], tleap_in_path])
        if tleap_out_path is None:
            tleap_out_path = fs.get_valid_temp_name(f"{eh_config['system.SCRATCH_DIR']}/tleap.out")
            temp_path_list.append(tleap_out_path)
        # write tleap.in
        with open(tleap_in_path, "w") as of:
            of.write(tleap_in_str)
        # run tleap command
        cmd_args = f"-f {tleap_in_path} > {tleap_out_path}"
        if if_ignore_start_up:
            cmd_args = f"-s {cmd_args}"
        if additional_search_path:
            for add_path in additional_search_path:
                cmd_args = f"{cmd_args} -I {add_path}"
        try:
            self.env_manager_.run_command("tleap", cmd_args)
        except CalledProcessError as e:
            if (not e.stderr.strip()) and (not e.stdout.strip()):  # empty stderr & stdout
                # find the error information in tleap.out
                new_e = self._find_tleap_error(tleap_out_path)
                # express the error
                for e_info in new_e.error_info_list:
                    _LOGGER.error(e_info)
                raise new_e from e

        # clean up temp file if success
        fs.clean_temp_file_n_dir(temp_path_list)

    @staticmethod
    def _find_tleap_error(tleap_out_path: str) -> tLEaPError:
        """an internal used function that find the text describing the error
        from a tleap output file
        Return a tLEaPError containing all the error information"""
        error_info_pattern = r"(.*(?:FATAL|Fatal).*\n(?:.*\n)+?)(?:(?:Exiting LEaP:)|/)"
        with open(tleap_out_path) as f:
            error_info_list = re.findall(error_info_pattern, f.read())
            error_info_list = [i.strip() for i in error_info_list]
            if not error_info_list:
                _LOGGER.warning("did not found any error in tleap out file. Here is the complete output:")
                _LOGGER.warning(f.read())

        return tLEaPError(error_info_list)

    # -- antechamber/parmchk2 --
    def run_antechamber(self, in_file: str, out_file: str, net_charge: int, spin: int,
                        charge_method: str, atom_type: str = "gaff", 
                        charge_file: str=None,
                        res_name: str=None,) -> None:
        """the python wrapper of running antechamber
        Args:
            in_file: input file path
            out_file: the output molecule description file path (e.g.: mol2, prepin, ac)
            atom_type: the target atom types for generation.
                = gaff : the default
                = gaff2: for GAFF, version 2
                = amber: for PARM94/99/99SB
                = bcc : for AM1-BCC
                = sybyl: for atom types used in sybyl
            charge_method: the method for determining the atomic charge of atoms from the input file
            charge_file: the charge file. (used when charge_method="rc")
            net_charge: net_charge of the molecule. Always explicitly specify it.
            spin: multiplicity (2S+1). Always explicitly specify it.
            res_name: residue name (only used if not available in the input file)

        Unsupported options: (Add when need)
            -a additional file name
            -fa additional file format
            -ao additional file operation
                crd : only read in coordinate
                crg: only read in charge
                radius: only read in radius
                name : only read in atom name
                type : only read in atom type
                bond : only read in bond type
            -rf residue topology file name in prep input file, default is molecule.res
            -ch check file name in gaussian input file, default is molecule
            -ek QM program (mopac or sqm) keyword (in quotes); overwrites previous keywords.
            -gk gaussian keyword in a pair of quotation marks
            -gm gaussian assign memory, inside a pair of quotes, such as "%mem=1000MB"
            -gn gaussian assign number of processor, inside a pair of quotes, such as "%nproc=8"
            -gv add keyword to generate gesp file (for Gaussian 09 only) 1: yes; 0: no, the default
            -ge gaussian esp file generated by iop(6/50=1), default is g09.gesp
            -df use divcon flag, 0 - use mopac; 2 - use sqm (the default)
            -du check atom name duplications, can be yes(y) or no(n), default is yes
            -bk 4-character component Id, for ccif
            -j atom type and bond type prediction index, default is 4
                0 : no assignment
                1 : atom type
                2 : full bond types
                3 : part bond types
                4 : atom and full bond type
                5 : atom and part bond type
            -eq equalize atomic charge, default is 1 for '-c resp' and '-c bcc'
                0 : no equalization
                1 : by atomic paths
                2 : by atomic paths and geometry, such as E/Z configurations
            -s status information, can be 0 (brief), 1 (the default) and 2 (verbose)
            -pf remove the intermediate files: can be yes (y) and no (n, default)
            -pl maximum path length to determin equivalence of atomic charges for resp and bcc.
                The smaller the value, the faster the algorithm, default is -1 (use full length),
                set this parameter to 10 to 30 if your molecule is big (# atoms >= 100)
                -dr acdoctor mode: validate the input file a la acdoctor, yes(y, default) or no(n)"""
        in_format = amber_interface.get_file_format(in_file)
        cmd_args = ["-i", in_file,
                    "-fi", in_format,
                    "-o", out_file,
                    "-fo", amber_interface.get_file_format(out_file),
                    "-nc", str(net_charge),
                    "-m", str(spin),
                    "-at", atom_type.lower()]
        if res_name:
            cmd_args.extend(["-rn", res_name])
        # charge
        if charge_method in self.SUPPORTED_CHARGE_METHOD_MAPPER:
            charge_method = self.SUPPORTED_CHARGE_METHOD_MAPPER[charge_method]
            cmd_args.extend(["-c", charge_method])
            if charge_method == "rc":
                if charge_file:
                    cmd_args.extend(["-cf", charge_file])
                else:
                    _LOGGER.error("charge_method='rc' is used but no charge_file is provided.")
                    raise ValueError
            if charge_method == "resp":
                if in_format not in ["gesp", "gout"]:
                    _LOGGER.error(f"charge_method='resp' is used but input file format ({in_file}) is not gesp or gout.")
                    raise ValueError
        else:
            _LOGGER.error(f"found unsupported charge method {charge_method}."
                          f"(Supported keywords: {self.SUPPORTED_CHARGE_METHOD_MAPPER.keys()})"
                          "Contact author for support if it is necessary for your work.")
            raise ValueError

        self.env_manager_.run_command("antechamber", cmd_args)

        # clean up
        fs.clean_temp_file_n_dir([
            "ATOMTYPE.INF",
            "NEWPDB.PDB",
            "PREP.INF",
            "sqm.pdb",
            "sqm.in",
            "sqm.out",
        ] + glob.glob("ANTECHAMBER*"))

    def run_parmchk2(self, in_file: str, out_file: str, gaff_type: str,
                     custom_force_field: str=None,
                     is_custom_ff_type_amber: bool=True,
                     ) -> None:
        """the python wrapper of running parmchk2
        Args:
            in_file: input file path
            out_file: the output frcmod file path
            gaff_type: ff parm set, it is suppressed by "custom_force_field" option
            custom_ff: the path of a customize ff file for search.
                       (e.g.: this allows you to use ff14SB for maa)
            is_custom_ff_type_amber: is the type of custom_ff amber
            * TODO: support -c -a -w when needed"""
        cmd_args = ["-i", in_file,
                    "-f", amber_interface.get_file_format(in_file),
                    "-o", out_file,
                    "-s", gaff_type.lower(),]
        if custom_force_field:
            cmd_args.extend(["-p", custom_force_field])
            if not is_custom_ff_type_amber:
                cmd_args.extend(["-pf", "2"])

        self.env_manager_.run_command("parmchk2", cmd_args)

    # -- cpptraj --

    # -- sander --

    # -- pmemd --

    # == engines ==
    # (engines for science APIs)
    # -- mutation/clean up --
    def tleap_clean_up_stru(
        self,
        input_pdb: str,
        out_path: str,
        if_align_index: bool = True,
        amber_lib: str = "leaprc.protein.ff14SB",
    ) -> None:
        """Method that uses tleap to clean up {input_pdb} by loading and saving the PDB with
        the {amber_lib}.
        Typical changes are:
            - complete missing atoms based on the AA defination in {amber_lib}
            of the correponding residue name
            - change terminal atom names
            - (if_align_index = False)renumber residue index and atom index from 1 and
            ignore chain index
            - add TER around all residue unit that is not recognized or is a ligand.
            (The might be harmful when a structure have a modified animo acid TODO)

        Args:
            input_pdb:
                the path of the input pdb for cleaning
            out_path:
                the path of the output pdb after cleaning
            if_align_index:
                whether or not align index back after cleaning.
            amber_lib:
                the amber library used for cleaning.

        Application:
            Used in mutation.api.mutate_stru_with_tleap()"""

        tleap_in_lines: List[str] = [
            f"source {amber_lib}",
            f"a = loadpdb {input_pdb}",
            f"savepdb a {out_path}",
            "quit",
        ]
        tleap_in_str = "\n".join(tleap_in_lines)
        self.run_tleap(tleap_in_str)
        if if_align_index:
            # temp file
            temp_dir = eh_config["system.SCRATCH_DIR"]
            fs.safe_mkdir(temp_dir)
            renumbered_pdb = fs.get_valid_temp_name(f"{temp_dir}/tleap_out_renumbered.pdb")
            # renumber
            index_mapper = pdb_io.get_index_mapper_from_pdb(out_path, input_pdb, method="by_order")
            pdb_io.restore_pdb_index(index_mapper, out_path, renumbered_pdb)
            shutil.move(renumbered_pdb, out_path)
            fs.clean_temp_file_n_dir([temp_dir, renumbered_pdb])

    # -- parameterize --
    def build_md_parameterizer(
            self,
            # find default values in AmberConfig
            # method
            force_fields: List[str] = "default",
            charge_method: str = "default",
            resp_engine: str = "default",
            resp_lvl_of_theory: str = "default",
            # workflow
            ncaa_param_lib_path: str = "default",
            force_renew_ncaa_parameter: bool = "default",
            ncaa_net_charge_engine: str = "default",
            ncaa_net_charge_ph: float = "default",
            solvate_box_type: str = "default",
            solvate_box_size: float = "default",
            # mod_aa_engine TODO ask if the other way really work
            # MCPY.py related TODO add when use
            # GB related
            gb_radii: str = None,
            # temp paths
            parameterizer_temp_dir: str = "default", # default: {config[system.SCRATCH_DIR]}/amber_parameterizer
            additional_tleap_lines: List[str] = None,
            ) -> AmberParameterizer:
        """the constructor for AmberParameterizer
        Args:
            force_fields:
                The list of force fields used for parameterization in Amber tleap format
                (e.g.: ["leaprc.protein.ff14SB", "leaprc.gaff", "leaprc.water.tip3p"])
            charge_method:
                The method used for determine the atomic charge.
                This method is applied to parameterization of ligand, modified AA,
                and metal binding site.
            (the following args only apply when charge_method="RESP")
                resp_engine:
                    The engine for calculating the RESP charge.
                resp_lvl_of_theory:
                    The level of theory for calculating the RESP charge
            ncaa_param_lib_path:
                The path of the non-CAA parameter library. This is where all generated
                NCAA params goes to. It will prevent redundant generation of same NCAAs.
                Normally we suggest setting this to a directory that contains all workflows
                of a same wild-type/template enzyme.
                * The NCAA-file correspondence is determined by the 
                  (1) 3-letter name in the file
                  (2) (if not 1 not exist) the file name
                * Setting this to a path that is too general may cause conflict when different
                NCAAs have the same name. (e.g.: different tautomer or general res name like LIG)
            force_renew_ncaa_parameter:
                Whether force renew the parameter files (frcmod etc.) for all ncaa
                (ligand, mod AA, or metal)
            ncaa_net_charge_engine:
                The engine the determines the net charge of NCAA if none is assigned in NCAA objects (Ligand,
                ModifedResidue TODO, MetalUnit)
            ncaa_net_charge_ph:
                The pH value used in determining the net charge of NCAA.
            solvate_box_type:
                The shape of the solvation box.
            solvate_box_size:
                The size of the solvation box.
            gb_radii:
                The effective GB radii used in the Generalized Born calculation. This will influence
                the GB radii in the prmtop file and are only used implicit solvent calculations.
            parameterizer_temp_dir: (default: {config[system.SCRATCH_DIR]}/amber_parameterizer)
                The temporary working directory that contains all the files generated by the AmberParameterizer
            additional_tleap_lines:
                handle for adding additional tleap lines before generating the parameters."""
        # tool: write the below code
        # print(AmberInterface._generate_default_assigning_lines_for_build_md_parameterizer(locals().items()))

        # init default values
        type_hint_sticker: AmberConfig
        if force_fields == "default":
            force_fields = self.config()["DEFAULT_FORCE_FIELDS"]        
        if charge_method == "default":
            charge_method = self.config()["DEFAULT_CHARGE_METHOD"]
        if resp_engine == "default":
            resp_engine = self.config()["DEFAULT_RESP_ENGINE"]
        if resp_lvl_of_theory == "default":
            resp_lvl_of_theory = self.config()["DEFAULT_RESP_LVL_OF_THEORY"]
        if ncaa_param_lib_path == "default":
            ncaa_param_lib_path = self.config()["DEFAULT_NCAA_PARAM_LIB_PATH"]
        if force_renew_ncaa_parameter == "default":
            force_renew_ncaa_parameter = self.config()["DEFAULT_FORCE_RENEW_NCAA_PARAMETER"]
        if ncaa_net_charge_engine == "default":
            ncaa_net_charge_engine = self.config()["DEFAULT_NCAA_NET_CHARGE_ENGINE"]
        if ncaa_net_charge_ph == "default":
            ncaa_net_charge_ph = self.config()["DEFAULT_NCAA_NET_CHARGE_PH"]
        if solvate_box_type == "default":
            solvate_box_type = self.config()["DEFAULT_SOLVATE_BOX_TYPE"]
        if solvate_box_size == "default":
            solvate_box_size = self.config()["DEFAULT_SOLVATE_BOX_SIZE"]
        if parameterizer_temp_dir == "default":
            parameterizer_temp_dir = self.config()["DEFAULT_PARAMETERIZER_TEMP_DIR"]

        return AmberParameterizer(
            force_fields,
            charge_method,
            resp_engine,
            resp_lvl_of_theory,
            ncaa_param_lib_path,
            force_renew_ncaa_parameter,
            ncaa_net_charge_engine,
            ncaa_net_charge_ph,
            solvate_box_type,
            solvate_box_size,
            gb_radii,
            parameterizer_temp_dir,
            additional_tleap_lines,)

    @staticmethod
    def _generate_default_assigning_lines_for_build_md_parameterizer(parameter_mapper: dict) -> str:
        """function only used for generate the code of the init default value section
        of build_md_parameterizer"""
        result = ""
        for k, v in parameter_mapper:
            if v == "default":
                result += (f"""if {k} == "default":
            {k} = self.config()["DEFAULT_{k.upper()}"]
        """)
        return result

    def antechamber_ncaa_to_moldesc(self,
                                    ncaa: NonCanonicalBase,
                                    out_path: str,
                                    gaff_type: str = "GAFF",
                                    charge_method: str = "AM1BCC",
                                    cluster_job_config: Dict=None,) -> str:
        """use antechamber to generate .mol2/.ac file for ligand/modified amino acid.
        Args:
            ncaa: the target Ligand/ModifiedAminoAcid
            out_path: the path of the output molecule description file. (use ext here to determine target format)
            gaff_type: the ff type used for NCAA. This influence the atom type in the moldesc file
            charge_method: the method for generation of atomic charges
            cluster_job_config: the configuration for submitting cluster jobs (used for Gaussian if charge_method='resp')
        Return:
            the out_path
            """
        # san check
        if (not ncaa.multiplicity) or (not ncaa.net_charge):
            _LOGGER.error("supplied NCAA does not have charge and spin."
                          " ALWAYS check and explicit assign it using"
                          " Structure.assign_ncaa_chargespin()")
            raise ValueError

        # 1. make ligand PDB
        temp_dir = eh_config["system.SCRATCH_DIR"]
        fs.safe_mkdir(temp_dir)
        temp_pdb_path = fs.get_valid_temp_name(f"{temp_dir}/{ncaa.name}.pdb")
        pdb_io.PDBParser().save_structure(temp_pdb_path, ncaa)
        input_file = temp_pdb_path

        # 2. determine charge
        charge_method = self.SUPPORTED_CHARGE_METHOD_MAPPER.get(charge_method, None)
        if charge_method:
            if charge_method == "bcc":
                pass
            if charge_method == "resp":
            # 2.1. Run RESP calculation (option)
                gout_file = self._calculate_charge_with(charge_method, ncaa, cluster_job_config)
                input_file = gout_file
        else:
            _LOGGER.error(f"found unsupported charge method {charge_method}."
                          f"(Supported keywords: {self.SUPPORTED_CHARGE_METHOD_MAPPER.keys()})")
            raise ValueError

        # 3. run antechamber on the PDB get mol2
        self.run_antechamber(in_file=input_file,
                             out_file=out_path,
                             charge_method=charge_method,
                             spin=ncaa.multiplicity,
                             net_charge=ncaa.net_charge,
                             atom_type=gaff_type,)

        # 4. clean up
        fs.clean_temp_file_n_dir([
            input_file,
            temp_dir
        ])
        return out_path

    def _calculate_charge_with(self, charge_method: str, ncaa: NonCanonicalBase,
                               cluster_job_config: Dict) -> str:
        """generate the charge file of {charge_method} for antechamber_ncaa_to_moldesc.
        return the output file containing the charge information varied by charge_method."""
        temp_chg_file = fs.get_valid_temp_name(
            f"{eh_config['system.SCRATCH_DIR']}/temp_charge.out")
        temp_opt_file = fs.get_valid_temp_name(
            f"{eh_config['system.SCRATCH_DIR']}/temp_opt.out")
        supported_list = ["resp"]
        if charge_method in supported_list:
            if charge_method == "resp":
                # 1. optimize (H only) "PBE1PBE/def2SVP em=gd3"? Do we?
                gaussian_interface.gaussain_optimize(ncaa,
                    out_file=temp_opt_file,
                    method="PBE1PBE/def2SVP em=gd3 nosymm",
                    cluster_job_config=cluster_job_config)
                # 2. calculate charge
                gaussian_interface.gaussain_single_point(temp_opt_file,
                    out_file=temp_chg_file,
                    method="hf/6-31g* SCF=tight Pop=MK iop(6/33=2)",
                    addition_output=True,
                    cluster_job_config=cluster_job_config,)
                # 3. clean up
                fs.clean_temp_file_n_dir(temp_opt_file)
        else:
            _LOGGER.error(f"using unsupported charge_method {charge_method}.")
            raise ValueError

        return temp_chg_file


    # region == TODO ==
    def write_minimize_input_file(self, fname: str, cycle: int) -> None:
        """Creates a minimization file to be used in an amber run. SHOULD NOT BE CALLED BY USERS DIRECTLY.
        All parameters in the &ctrl block are hardcoded except for ncyc and ntpr, which are 0.5*cycle
        and 0.2*cycle as integers, respectively.
        """
        minimize_lines: List[str] = [
            "Minimize",
            " &cntrl",
            "  imin=1,",
            "  ntx=1,",
            "  irest=0,",
            f"  maxcyc={str(cycle)},",
            f"  ncyc={int(0.5 * cycle)},",
            f"  ntpr={int(0.2 * cycle)},",
            "  ntwx=0,",
            "  cut=8.0,",
            " /",
            "",
        ]
        fs.write_lines(fname, minimize_lines)

    def minimize_structure(self, pdb: str, min_dir: str = "./", mode: str = "CPU", cycle: int = 2000) -> str:
        """Class method that minimizes the structure found in a supplied .pdb file, returning
        the path to the minimized file.
        Args:
            pdb: The .pdb file with the structure to be minimized.
            mode: Which version of Amber to use. Allowed values are "CPU" and "GPU".
            min_dir: Work directory where the paramter files and output are saved.
            cycle: Number of steps in the steepest descent of the cycle.
        Returns:
            Path to the minimized structure in a .pdb file.
        """
        # TODO(CJ): add some checks for inputs
        inpath = Path(pdb)
        if min_dir == "./":
            min_dir = f"{inpath.parent}/pdb_min/"

        fs.safe_mkdir(min_dir)
        outfile = f"{min_dir}/{inpath.stem}_min.pdb"
        min_in = f"{min_dir}/min.in"
        min_out = f"{min_dir}/min.out"
        min_rst = f"{min_dir}/min.ncrst"
        (prmtop, inpcrd) = self.build_param_files(str(inpath), min_dir)
        self.write_minimize_input_file(min_in, cycle)
        engine = self.config_.get_engine(mode)

        # yapf: disable
        self.env_manager_.run_command(
            engine, ["-O", "-i", min_in, "-o", min_out, "-p", prmtop, "-c", inpcrd, "-r", min_rst],

        )

        self.env_manager_.run_command(
            "ambpdb",[f"-p", prmtop, "-c", min_rst, ">", outfile]
        )
        # yapf: enable

        return outfile

    def _setup_solvation(self, igb: Union[None, str]) -> List[str]:
        """Method that creates solvation box settings for system paramterization.
        Args:
            igb: The Born solvation setting.
        Returns:
            A list() of the lines to be added to the leap.in. 
        """
        result: List[str] = list()
        if igb:
            igb_str: str = str(igb)
            radii: str = self.config_.RADII_MAP.get(igb_val)
            if radii:
                result.append(f"set default PDBRadii {radii}")
            else:
                _LOGGGER.warning(f"The IGB value {igb} is NOT suppported. Continuing...")

        result.append("center a")
        result.extend(["addions a Na+ 0", "addions a Cl- 0"])

        if self.config_.BOX_TYPE == "oct":
            result.append(f"solvateOct a TIP3PBOX {self.config_.BOX_SIZE}")
        elif self.config_.BOX_TYPE == "box":
            result.append(f"solvatebox a TIP3PBOX {self.config_.BOX_SIZE}")
        else:
            _LOGGER.error(f"The supplied box type {self.config_.BOX_SIZE} is NOT supported. Exiting...")
            exit(1)
        return result

    def build_param_files(self, in_pdb: str, build_dir: str, pH: float = 7.0, igb: str = None) -> Tuple[str, str]:
        """Creates the .prmtop and .inpcrd files for the supplied .pdb file. Handles
        processing of the Ligand() and MetalCenter() objects in the structure. Note that the supplied
        pdb file is assumed to be protonated.
        Args:
            in_pdb: The .pdb file to build parameter files for.
            buld_dir: The directory to build the parameter files in.
            pH: A float() representing the pH the parameter files should be made at. Default value is 7.0.
            igb: a str() representing Born sovlation radius settting. 
        Returns:
            A Tuple[str,str] with the containing (.prmtop path, .inpcrd path).
        """
        ligand_dir: str = f"{build_dir}/ligands/"
        metalcenter_dir: str = f"{build_dir}/metalcenter/"

        fs.safe_mkdir(ligand_dir)
        fs.safe_mkdir(metalcenter_dir)

        structure: Structure = PDBParser().get_structure(in_pdb)
        ligand_paths: List[str] = structure.build_ligands(ligand_dir, True)
        for lig in ligand_paths:
            _ = prep.protonate._protonate_ligand_PYBEL(lig, pH, lig)
        ligand_charges: List[int] = list(map(lambda pp: prep.protonate._ob_pdb_charge(pp), ligand_paths))
        ligand_params: List[Tuple[str, str]] = self.build_ligand_param_files(ligand_paths, ligand_charges)
        leap_path: str = f"{build_dir}/leap.in"
        leap_log: str = f"{build_dir}/leap.out"

        leap_contents: List[str] = [
            "source leaprc.protein.ff14SB",
            "source leaprc.gaff",
            "source leaprc.water.tip3p",
        ]
        for (prepin, frcmod) in ligand_params:
            leap_contents.extend([f"loadAmberPrep {prepin}", f"loadAmberParams {frcmod}"])
        leap_contents.append(f"a = loadpdb {in_pdb}")

        leap_contents.extend(self._setup_solvation(igb))
        pdb_path: Path = Path(in_pdb)
        prmtop: str = f"{build_dir}/{pdb_path.stem}.prmtop"
        inpcrd: str = f"{build_dir}/{pdb_path.stem}.inpcrd"
        pdb_ff: str = f"{build_dir}/{pdb_path.stem}_ff.pdb"
        leap_contents.extend([f"saveamberparm a {prmtop} {inpcrd}", f"savepdb a {pdb_ff}", "quit"])
        fs.write_lines(leap_path, leap_contents)

        # yapf: disable
        self.env_manager_.run_command(
            "tleap",
            ["-s", "-f", leap_path, ">", leap_log]
        )
        # yapf: enable

        good_output: bool = True

        for pp in map(Path, [prmtop, inpcrd, pdb_ff]):
            msg = str()
            if not pp.exists():
                good_output = False
                msg = f"The file {pp} does not exist."
            else:
                if pp.stat().st_size == 0:
                    msg = f"The file {pp} exists but is empty."
                    good_output = False
            if msg:
                _LOGGER.warn(msg)

        if not good_output:
            _LOGGER.error(f"The output generated by AmberInterface.build_param_files() is incomplete. Exiting...")
            exit(1)

        return (prmtop, inpcrd)

    def build_ligand_param_files(
        self,
        paths: List[str],
        charges: List[int],
    ) -> List[Tuple[str, str]]:
        # TODO(CJ): add the method flag?
        """Creates .prepin and .frcmod files for all the supplied .pdb files. Saves files to
        same directory as the supplied .pdb files. Removes intermediate files. Should not
        be called directly by the user. Instead use AmberInterface.build_param_files()
        Args:
            paths: A list() of ligand .pdb files to prepare. MUST BE SAME LENGHT AS charges.
            charges: A list() of charges accompanying the paths. MUST BE SAME LENGTH AS paths.
        Returns:
            A list() of filename pairs with with the .prepin and .frcmod files for each supplied .pdb file.
        Raises:
            AssertionErrors: When various input sanitization checks fail.
        """
        result: List[Tuple[str, str]] = list()
        assert len(paths) == len(charges)
        for lig_pdb, net_charge in zip(paths, charges):
            lig_name: str = PDBParser.get_structure(lig_pdb)
            print(lig_name)
            lig_pdb = Path(lig_pdb)
            prepin: str = str(lig_pdb.with_suffix(".prepin"))
            frcmod: str = str(lig_pdb.with_suffix(".frcmod"))
            lig_pdb: str = str(lig_pdb)
            # TODO(CJ): check if you can get the ligand name from the
            # .pdb filename alone... I think this may be possible
            # if renew
            # yapf: disable
            self.env_manager_.run_command(
                "antechamber",
                ["-i", lig_pdb, "-fi", "pdb", "-o", prepin, "-fo", "prepi", "-c", "bcc", "-s", "0", "-nc", str(net_charge), ],
            )
            # yapf: enable
            self.remove_antechamber_temp_files()
            self.env_manager_.run_command(
                # yapf: disable
                "parmchk2",
                ["-i", prepin, "-f", "prepi", "-o", frcmod])
            # yapf: enable
            # record
            result.append((prepin, frcmod))
        return result

    def md_min_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_MIN, creates a min.in file for an Amber minimization run.
        These settings are updated by accessing the owned AmberConfig object through AmberInterface.config().
        Args:
            outfile: The name of the file to save the minimization input file to..
        Returns:
            Path to the minimization file.
        """
        config = self.config_.CONF_MIN

        contents: List[str] = [
            "Minimize",
            " &cntrl",
            "  imin  = 1,",
            "  ntx   = 1,",
            "  irest = 0,",
            f"  ntc = {config['ntc']},",
            f"  ntf = {config['ntf']},",
            f"  cut   = {config['cut']},",
            f"  maxcyc = {config['maxcyc']},",
            f"  ncyc  = {int(config['ncyc_mult']*config['maxcyc'])},",
            f"  ntpr = {int(config['ntpr_mult']*config['maxcyc'])},",
            f"  ntwx = 0,",
        ]

        last_line: str = ""
        if config["ntr"] == 1:
            last_line = f"  ntr   = {config['ntr']}, restraint_wt = {config['restraint_wt']}, restraintmask = {config['restraintmask']},"

        contents.append(last_line)
        contents.extend(["/", ""])
        fs.write_lines(outfile, contents)
        return outfile

    def md_heat_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_HEAT, creates a heat.in file for an Amber heating run.
        These settings are updated by accessing the owned AmberConfig object through AmberInterface.config().
        Args:
            outfile: The name of the file to save the heating input file to.
        Returns:
            Path to the heating input file.
        """
        config = self.config_.CONF_HEAT

        contents: List[str] = [
            "Heat",
            " &cntrl",
            "  imin = 0,",
            "  ntx = 1,",
            "  irest = 0,",
            f"  ntc = {config['ntc']},",
            f"  ntf = {config['ntf']},",
            f"  cut = {config['cut']},",
            f"  nstlim = {config['nstlim']},",
            f"  dt = {config['dt']},",
            f"  tempi = {config['tempi']},",
            f"  temp0 = {config['temp0']},",
            f"  ntpr =  {int(config['ntpr']*config['nstlim'])},",
            f"  ntwx =  {int(config['ntwx']*config['nstlim'])},",
            f"  ntt = {config['ntt']},",
            f"  gamma_ln = {config['gamma_ln']},",
            "  ntb = 1,",
            "  ntp = 0,",
            f"  iwrap = {config['iwrap']},",
            "  nmropt = 1,",
            "  ig = -1,",
        ]

        if config["ntr"] == 1:
            contents.append(f"  ntr = {config['ntr']},")
            contents.append(f"  restraint_wt = {config['restraint_wt']},")
            contents.append(f"  restraintmask = {config['restraintmask']},")
        contents[-1] = contents[-1][0:-1]
        contents.extend([
            " /",
            " &wt",
            "  type = 'TEMP0',",
            "  istep1 = 0,",
            f"  istep2 = {int(config['A_istep2']*config['nstlim'])},",
            f"  value1 = {config['tempi']},",
            f"  value2 = {config['temp0']},",
            " /",
            f" &wt",
            "  type = 'TEMP0',",
            f"  istep1 = {int(config['A_istep2']*config['nstlim'])+1},",
            f"  istep2 = {config['nstlim']},",
            f"  value1 = {config['temp0']},",
            f"  value2 = {config['temp0']},",
            " /",
            " &wt",
            "  type = 'END',",
            " /",
            "",
        ])

        fs.write_lines(outfile, contents)
        return outfile

    def md_equi_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_EQUI, creates an input file for an Amber constant
        equilibration run. These settings are updated by accessing the owned AmberConfig object through
        AmberInterface.config().
        Args:
            outfile: The name of the file to save the constant pressure equilibration input.
        Returns:
            Path to the constant pressure equilibration input file.
        """
        config = self.config_.CONF_EQUI
        contents: List[str] = [
            "Equilibration:constant pressure",
            " &cntrl",
            "  imin = 0,",
            f"  ntx = {config['ntx']},",
            f"  irest = {config['irest']},",
            f"  ntf = {config['ntf']},",
            f"  ntc = {config['ntc']},",
            f"  nstlim = {config['nstlim']},",
            f"  dt = {config['dt']},",
            f"  cut = {config['cut']},",
            f"  temp0 = {config['temp0']},",
            f"  ntpr = {int(config['ntpr']*config['nstlim'])},",
            f"  ntwx = {config['ntwx']},",
            f"  ntt = {config['ntt']},",
            f"  gamma_ln = {config['gamma_ln']},",
            "  ntb = 2,",
            "  ntp = 1,",
            f"  iwrap = {config['iwrap']},",
            "  ig = -1,",
        ]

        if config["ntr"] == 1:
            contents.extend([
                f"  ntr   = {config['ntr']},",
                f"  restraint_wt = {config['restraint_wt']},",
                f"  restraintmask = {config['restraintmask']},",
            ])

        contents.extend(["/", ""])

        fs.write_lines(outfile, contents)
        return outfile

    def md_prod_file(self, outfile: str) -> str:
        """Using the settings specified by AmberConfig.CONF_PROD, creates an input file for an Amber production
        md run. These settings are updated by accessing the owned AmberConfig object through
        AmberInterface.config().
        Args:
            outfile: The name of the file to save the production md input.
        Returns:
            Path to the production md input file.
        """
        config = self.config_.CONF_PROD
        contents: List[str] = [
            "Production: constant pressure",
            " &cntrl",
            "  imin = 0,",
            f"  ntx = {config['ntx']},",
            f"  irest = {config['irest']},",
            f"  ntf = {config['ntf']},",
            f"  ntc = {config['ntc']},",
            f"  nstlim = {config['nstlim']},",
            f"  dt = {config['dt']},",
            f"  cut = {config['cut']},",
            f"  ntpr = {int(config['ntpr']*config['nstlim'])},",
            f"  ntwx = {config['ntwx']},",
            f"  ntt = {config['ntt']},",
            f"  gamma_ln = {config['gamma_ln']},",
            "  ntb = 2,",
            "  ntp = 1,",
            f"  iwrap = {config['iwrap']},",
            "  ig = -1,",
        ]

        if config["ntr"] == 1:
            contents.extend([
                f"  ntr = {config['ntr']},",
                f"  restraint_wt = {config['restraint_wt']},",
                f"  restraintmask = {config['restraintmask']}",
            ])

        contents.extend(["/", ""])
        fs.write_lines(outfile, contents)
        return outfile

    def md_run(self, prmtop: str, inpcrd: str, work_dir: str, mode: str = "CPU") -> str:
        """Runs a full MD simulation using the suplied prmtop and inpcrd files. Simulation is composed
        of four steps:
                1. minimization
                2. heating
                3. equilibration
                4. production
        The parameter files for these steps are created by AmberInterface.md_min_file(), AmberInterface.md_heat_file(),
        AmberInterface.md_equi_file() and AmberInterface.md_prod_file(), respectively. All work is done in the supplied
        work_dir location.
        Args:
            prmtop: str() with a path to a .prmtop file generated by AmberInterface.build_param_files().
            inpcrd: str() with a path to a .inpcrd file generated by AmberInterface.build_param_files().
            work_dir: Directory where the temporary files should be stored.
            mode: The mode to use for calculations. Only "CPU" and "GPU" are allowed.
        Returns:
            Path to the production .nc file.
        """
        fs.safe_mkdir(work_dir)

        min_in: str = self.md_min_file(f"{work_dir}/min.in")
        heat_in: str = self.md_heat_file(f"{work_dir}/heat.in")
        equi_in: str = self.md_equi_file(f"{work_dir}/equi.in")
        prod_in: str = self.md_prod_file(f"{work_dir}/prod.in")
        engine = self.config_.get_engine(mode)

        min_out: str = f"{work_dir}/min.out"
        min_rst: str = f"{work_dir}/min.rst"

        #yapf: disable
        self.env_manager_.run_command(
            engine,
            ["-O", "-i", min_in, "-o", min_out, "-p", prmtop, "-c", inpcrd, "-r", min_rst, "-ref", inpcrd, ],
        )
        #yapf: enable

        heat_out: str = f"{work_dir}/heat.out"
        heat_rst: str = f"{work_dir}/heat.rst"

        #yapf: disable
        self.env_manager_.run_command(
            engine,
            ["-O", "-i", heat_in, "-o", heat_out, "-p", prmtop, "-c", min_rst, "-ref", min_rst, "-r", heat_rst, ],
        )
        #yapf: enable

        equi_out: str = f"{work_dir}/equi.out"
        equi_rst: str = f"{work_dir}/equi.rst"
        equi_nc: str = f"{work_dir}/equi.nc"

        # yapf: disable
        self.env_manager_.run_command(
            engine,
            ["-O", "-i", equi_in, "-o", equi_out, "-p", prmtop, "-c", heat_rst, "-ref", heat_rst, "-r", equi_rst, "-x", equi_nc, ],
        )
        # yapf: enable

        prod_out: str = f"{work_dir}/prod.out"
        prod_rst: str = f"{work_dir}/prod.rst"
        prod_nc: str = f"{work_dir}/prod.nc"

        # yapf: disable
        self.env_manager_.run_command(
            engine,
            ["-O",  "-i",  prod_in,  "-o",  prod_out,  "-p",  prmtop,  "-c",  equi_rst,  "-ref",  equi_rst,  "-r",  prod_rst,  "-x",  prod_nc, ],
        )
        # yapf: enable

        return prod_nc


    def nc2mdcrd(
        self,
        nc_in: str,
        prmtop: str,
        point: Union[int, None] = None,
        start: int = 1,
        end: Union[int, str] = "last",
        step: int = 1,
        engine: str = "cpptraj",
    ) -> str:
        """Converts the supplied .nc file to an .mdcrd file. Scope of trajectory conversion can be
        specified though by default all available frames are selected.
        Args:
            nc_in: The path to the .nc file as a str().
            prmtop: The path to the prmtop file as  str().
            point: If not None, step is (nstlim/ntwx)/point. Default is None, but value should be int() when given.
            start: 1-indexed starting point. Default value is 1.
            end: 1-indexed ending point. By default uses "last"
            engine: The engine to convert the .nc file. Default value is "cpptraj".
        Returns:
            The name of the converted .mdcrd file.
        """
        mdcrd: str = str(Path(nc_in).with_suffix(".mdcrd"))
        config = self.config_.CONF_PROD
        if point is not None:
            step: int = int((int(config["nstlim"]) / int(config["ntwx"])) / point)

        if engine == "cpptraj":
            cpptraj_in = "./cpptraj_nc2mdcrd.in"
            cpptraj_out = "./cpptraj_nc2mdcrd.out"
            contents: List[str] = [
                f"parm {prmtop}",
                f"trajin {nc_in} {start} {end} {step}",
                f"trajout {mdcrd}",
                "run",
                "quit",
            ]
            fs.write_lines(cpptraj_in, contents)
            self.env_manager_.run_command("cpptraj", ["-i", cpptraj_in, ">", cpptraj_out])
            #fs.safe_rm(cpptraj_in)
            #fs.safe_rm(cpptraj_out)
        else:
            raise UnsupportedMethod(f"'{engine}' is not a supported method for AmberInteface.nc2mdcrd(). Allowed methods are: 'cpptraj'.")

        return mdcrd

    def get_frames(
        self,
        nc_in: str,
        prmtop: str,
        point: Union[int, None] = None,
        start: int = 1,
        end: Union[int, str] = "last",
        step: int = 1,
    ) -> str:
        """Converts the supplied .nc file to an .mdcrd file. Scope of trajectory conversion can be
        TODO(CJ): update this
        """
        outfile: str = f"{Path(prmtop).parent}/cpptraj_frames.pdb"
        config = self.config_.CONF_PROD
        if point is not None:
            step: int = int((int(config["nstlim"]) / int(config["ntwx"])) / point)

        cpptraj_in = "./cpptraj_nc2mdcrd.in"
        cpptraj_out = "./cpptraj_nc2mdcrd.out"
        contents: List[str] = [
            f"parm {prmtop}",
            f"trajin {nc_in} {start} {end} {step}",
            "strip :WAT",
            "strip :Na+,Cl-",
            f"trajout {outfile} conect",
            "run",
            "quit",
        ]
        fs.write_lines(cpptraj_in, contents)
        self.env_manager_.run_command("cpptraj", ["-i", cpptraj_in, ">", cpptraj_out])
        #fs.safe_rm(cpptraj_in)
        fs.safe_rm(cpptraj_out)
        charges = read_charge_list(prmtop)
        result = frames_from_pdb(outfile)

        if not len(result):
            return result

        num_atoms = len(result[0].atoms)
        charges = charges[0:num_atoms]
        for fidx, frame in enumerate(result):
            result[fidx].update_charges(charges)

        return result

    def parse_fmt(self, fmt: str) -> Tuple[Any, int]:
        """Private helper function that parses a %FORMAT() string from an amber .prmtop file.
        Identifies both the type and width of the format. Supplied string can contain whitespace.
        In the case that the supplied format string does not contain the pattern '%FORMAT(<>)',
        the result (None, -1) is returned.

        Args:
            fmt: the %FORMAT() line from a .prmtop as a str().

        Returns:
            A tuple() of length 2 containing the (type ctor, width as an int()).
        """

        if fmt.find('%FORMAT(') == -1 or fmt.find(')') == -1:
            return (None, -1)

        fmt = fmt.strip().replace('%FORMAT(', '').replace(')', '')

        idx = 0

        while fmt[idx].isnumeric():
            idx += 1

        dtype = None

        type_string = fmt[idx].upper()

        if type_string == 'A':
            dtype = str
        elif type_string == 'I':
            dtype = int
        elif type_string == 'E':
            dtype = float

        width = int(fmt[idx + 1:].split('.')[0])

        return (dtype, width)

    def parse_prmtop(self, prmtop: str) -> Dict:
        """Parses a .prmtop file into a dict() with the pattern of (key, value) being (keyword after %FLAG, following values). 
        Takes the specified %FORMAT() into account. Checks if the supplied file exists, exiting if it does not.

        Args:
            prmtop: str() holding the path to a .prmtop file.

        Returns:
            A dict() with (key, value) pairs of (keyword, formatted/typed values).

        """
        #TODO(CJ): unit tests

        if not Path(prmtop).exists():
            _LOGGER.error(f"The file '{prmtop}' does not exist. Exiting...")
            exit(1)

        raw: str = fs.content_from_file(prmtop)

        result = dict()

        for chunk in raw.split('%FLAG'):
            if chunk.find('%VERSION') != -1:
                #TODO(CJ): figure out what I want to do with version info
                continue

            (key, fmt, body) = chunk.split('\n', 2)
            body = re.sub('\n', '', body)
            (dtype, width) = self.parse_fmt(fmt)

            raw: List[str] = list()
            for idx in range(0, len(body), width):
                raw.append(body[idx:idx + width].strip())

            result[key.strip()] = list(map(dtype, raw))

        return result

    def remove_antechamber_temp_files(self, dname: str = './') -> None:
        """Helper method that removes temporary files generated by antechamber from a given directory.
        Removes the files "ATOMTYPE.INF", "NEWPDB.PDB", "PREP.INF" and those with the patterns "ANTECHAMBER*"
        and "sqm.*".

        Args:
            dname: str() with the path to the directory to clean out. Uses current directory as default.

        Returns:
            Nothing.

        """

        files_to_remove: List[str] = [f"{dname}/ATOMTYPE.INF", f"{dname}/PREP.INF", f"{dname}/NEWPDB.PDB"]
        files_to_remove.extend(list(map(str, Path(f"{dname}").glob("ANTECHAMBER*"))))
        files_to_remove.extend(list(map(str, Path(f"{dname}").glob("sqm*"))))

        _ = list(map(fs.safe_rm, files_to_remove))

    def add_charges(self, stru: Structure, prmtop: str) -> None:
        """Method that adds RESP charges from a .prmtop file to the supplied Structure object. If the supplied prmtop
        file does not line up with the structure, an error is thrown. Performs operation in place.
        Will also throw an error if some of the resulting atoms do not get assigned charges.

        Args:
            stru: The Structure object that charges will be added to.
            prmtop: The str() path of an Amber paramter file.

        Returns:
            Nothing.
        """
        p_dict = self.parse_prmtop(prmtop)

        charges: List[float] = p_dict['CHARGE']
        anames: List[str] = p_dict['ATOM_NAME']
        res_pointers: List[int] = p_dict['RESIDUE_POINTER']
        res_names: List[int] = p_dict['RESIDUE_LABEL']

        temp = []
        for (ridx, rp), rn in zip(enumerate(res_pointers[:-1]), res_names):
            temp.extend([rn] * (res_pointers[ridx + 1] - rp))

        temp.extend((len(anames) - len(temp)) * [rn[-1]])
        res_names = temp

        for c, an, rn, aa in zip(charges, anames, res_names, stru.atoms):

            if not (an == aa.name and aa.residue.name == rn):
                _LOGGER.error(f"There is a mismatch in the atoms. Make sure the supplied prmtop file is correct. Exiting..")
                exit(1)
            aa.charge = c

        if not stru.has_charges():
            _LOGGER.error("Not all atoms in the supplied structure were assigned charges. Exiting...")
            exit(1)
            pass

    # endregion == TODO ==


amber_interface = AmberInterface(None, eh_config._amber)
