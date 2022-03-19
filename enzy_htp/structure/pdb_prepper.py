"""Defines a PDBPrepper class that takes in a raw pdb and prepares it for analysis.
For existing users of EnzyHTP, this largely replaces what Class_PDB did.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-03-12
"""
import shutil
import numpy as np

import pdb2pqr
from ..core.logger import _LOGGER, init_logger
from pdb2pqr.main import main_driver as run_pdb2pqr
from pdb2pqr.main import build_main_parser as build_pdb2pqr_parser

from typing import Set
from ..core.file_system import (
    base_file_name,
    get_current_time,
    safe_mkdir,
    lines_from_file,
    write_lines,
)
from .structure import structure_from_pdb

from .structure import Structure
from .pdb_line import PDBLine, read_pdb_lines

import openbabel
import openbabel.pybel as pybel

from .ligand import Ligand, protonate_ligand

from typing import List

# from .mutate import MutaFlag, mutaflag_to_str
from ..preparation import MutaFlag, mutaflag_to_str

# from ..core import em
from ..chemical import convert_to_three_letter, get_element_aliases


class PDBPrepper:
    def __init__(
        self, pdb_name, **kwargs
    ):  # , PDB_input, wk_dir="", name="", input_type="path"):
        """Inits PDBPrepper with a pdb file and optionally a work directory to place temporary files."""
        self.no_water_path = None
        self.pdb_path = pdb_name
        self.base_pdb_name = base_file_name(pdb_name)
        self.work_dir = kwargs.get(
            "work_dir", f"{get_current_time()}_{self.base_pdb_name}"
        )
        safe_mkdir(self.work_dir)
        shutil.copy(self.pdb_path, f"{self.work_dir}/{self.base_pdb_name}.pdb")
        self.path_name = f"{self.work_dir}/{self.base_pdb_name}.pdb"
        self.pqr_path = f"{self.work_dir}/{self.base_pdb_name}.pqr.pdb"
        self.mutations = []

    def rm_wat(self) -> str:
        """
        Remove water and ion for the pdb. Remians the same if there's no water or ion.
        Now only skip [Na+,Cl-,WAT,HOH] // Append more in the future.
        Save changed files into self.path.
        TODO: need to support key water.
        """

        pdb_lines: List[PDBLine] = read_pdb_lines(self.path_name)
        pdb_lines = list(
            filter(lambda pl: not (pl.is_water() or pl.is_CRYST1()), pdb_lines)
        )
        mask = [True] * len(pdb_lines)

        for pidx, pl in enumerate(pdb_lines[1:]):
            if pl.is_TER() and pdb_lines[pidx - 1].is_TER():
                mask[pidx] = False

        pdb_lines = np.array(pdb_lines)[mask]
        self.no_water_path = f"{self.work_dir}/{self.base_pdb_name}_rmW.pdb"
        write_lines(self.no_water_path, list(map(str, pdb_lines)))

        return self.no_water_path

    def _update_name(self):
        """
        update name
        """
        suffix_len = len(self.path.split(".")[-1]) + 1
        self.name = self.path.split(os.sep)[-1][:-suffix_len]
        self.path_name = self.dir + "/" + self.name

    def get_stru(self, ligand_list=None, renew=0):
        """
        Convert current PDB file (self.path) to a Struture object (self.stru).
        ------
        input_name  : a name tag for the object (self.name by default)
        ligand_list : a list of user assigned ligand residue names.
        renew       : 1: force generate a new stru_obj
        """
        # indicated by self.name
        input_name = self.name
        # if get new stru
        get_flag = 0
        if self.stru is not None:
            if self.stru.name != self.name:
                get_flag = 1
                # warn if possible wrong self.stru
                if Config.debug >= 1:
                    print("PDB.get_stru: WARNING: self.stru has a different name")
                    print("     -self.name: " + self.name)
                    print("     -self.stru.name: " + self.stru.name)
                    print("Getting new stru")
        else:
            get_flag = 1

        if get_flag or renew:
            if self.path is not None:
                self.stru = Structure.fromPDB(
                    self.path, input_name=input_name, ligand_list=ligand_list
                )
            else:
                self.stru = Structure.fromPDB(
                    self.file_str,
                    input_type="file_str",
                    input_name=input_name,
                    ligand_list=ligand_list,
                )

    def _get_file_str(self):
        """
        read file_str from path is not pre-exist.
        -------------
        recommend before every potential first use of self.file_str
        """
        if self.path is None:
            return self.file_str
        else:
            return open(self.path).read()

    def _get_file_path(self):
        """
        save a file and get path if self.path is None
        -------------
        recommend before every potential first use of self.path
        """
        if self.path is None:
            self.path = self.path_name + ".pdb"
            with open(self.path, "w") as of:
                of.write(self.file_str)
        return self.path

    def _init_MD_conf(self):
        """
        initialize default MD configuration. Replaced by manual setting if assigned later.
        """
        self.conf_min = Config.Amber.conf_min
        self.conf_heat = Config.Amber.conf_heat
        self.conf_equi = Config.Amber.conf_equi
        self.conf_prod = Config.Amber.conf_prod

    def set_oniom_layer(self, atom_list=[], preset=0):
        """
        set oniom layer with san check
        layer_atoms are in higher pirority
        """
        if len(atom_list) != 0:
            if all(type(x) is str for x in atom_list):
                self.layer_atoms = atom_list
            else:
                raise Exception(
                    "set_oniom_layer: atom_list require a list of str. e.g.: ['1-9','11,13-15']"
                )
        else:
            if preset == 0:
                raise Exception(
                    "set_oniom_layer: please assign one of the argument: atom_list (a list of layer_atoms) or preset (other than 0)"
                )
            else:
                self.layer_preset = preset

    def get_last_A_id(self):
        """
        get last atom id in PDB
        """
        with open(self.path) as f:
            fl = f.readlines()
            for i in range(len(fl) - 1, -1, -1):
                pdbl = PDB_line(fl[i])
                if pdbl.line_type == "ATOM" or pdbl.line_type == "HETATM":
                    return pdbl.atom_id

    """
    =========
    Sequence TODO: reform(most of it has been moved to class Chain，only judgement of art residue and overlook of whole structure remains)
    =========
    """

    def get_seq(self, Oneletter=0):
        """
        get_seq(self, Oneletter=0)
        Support most PDB types
        ----------------------------
        Oneletter
        = 0 // use 3-letter format to represet each residue
        = 1 // use 1-letter format to represet each residue
        ----------------------------
        Get sequence for current PDB (self.path) // A general function to obtain the missing residue
        + Use "NAN"/"-" as a filler to store missing residues (detect internal missing)
        + Check if contain non-standard residue with the Resi_map --> self.if_art_resi
        + Check if contain ligand with the Resi_map and TIP3P_map --> self.if_ligand
        - (WARNING) Require the ligand in seperate chains
        - Re-assign the chain_index base on the order in the file
        - Do not include the original residue index
        - Do not include any HETATM (ligand/solvent)
        ----------------------------
        save the info to self.sequence/self.sequence_one and return it
        ----------------------------
        self.sequence:
        Format: {'Chain_index':['res','res',...]
                 'Chain_index':['res','res',...]
                 ...
                }

        self.raw_sequence:
            Internal used, containing HETATM

        self.sequence_one: (when Oneletter=1)
        Format: {'Chain_index':'seq'
                 'Chain_index':'seq'
                 ...
                }

        !!NOTE!! the self.sequence needs to be updated after mutation
        """

        PDB_str = self._get_file_str()
        self.raw_sequence = {}
        self.sequence = {}
        self.sequence_one = {}

        Chain_str = PDB_str.split(line_feed + "TER")  # Note LF is required

        for chain, i in enumerate(Chain_str):

            Chain_index = chr(65 + chain)  # Covert to ABC using ACSII mapping
            Chain_sequence = []

            # Get the Chain_sequence
            lines = i.split(line_feed)

            for line in lines:

                pdb_l = PDB_line(line)

                if pdb_l.line_type == "ATOM" or pdb_l.line_type == "HETATM":

                    # Deal with the first residue
                    if len(Chain_sequence) == 0:
                        Chain_sequence.append(pdb_l.resi_name)
                        last_resi_index = pdb_l.resi_id
                        continue

                    # find next new residue
                    if pdb_l.resi_id != last_resi_index:

                        # Deal with missing residue, fill with "NAN"
                        missing_length = pdb_l.resi_id - last_resi_index - 1
                        if missing_length > 0:
                            Chain_sequence = Chain_sequence + ["NAN",] * missing_length

                        # Store the new resi
                        Chain_sequence.append(pdb_l.resi_name)

                    # Update for next loop
                    last_resi_index = pdb_l.resi_id

            self.raw_sequence[Chain_index] = Chain_sequence

        self._strip_raw_seq()  # strip the raw_sequence and save to sequence

        self.get_if_complete()

        if Oneletter == 1:
            self._get_Oneletter()
            return self.sequence_one
        else:
            return self.sequence

    get_if_ligand = get_seq  # alias
    get_if_art_resi = get_seq

    def _strip_raw_seq(self):
        """
        (Used internally) strip the raw_sequence.
        - Delete ligand and solvent
        - Delete chains without residue

        save changes to self.sequence

        Judge if containing any ligand or artificial residue

        save to self.if_ligand and self.if_art_resi
        """
        new_index = 0

        # if the value is not changed to 1 then it's 0
        self.if_art_resi = 0
        self.if_ligand = 0

        if len(self.raw_sequence) == 0:
            print("The self.raw_sequence should be obtained first")
            raise IndexError

        for chain in self.raw_sequence:

            chain_seq = []
            if_realchain = 0

            # Judge if real chain
            for name in self.raw_sequence[chain]:
                clean_name = name.strip(" ")
                if clean_name in Resi_map2.keys():
                    if_realchain = 1

            if if_realchain:
                for name in self.raw_sequence[chain]:
                    clean_name = name.strip(" ")
                    chain_seq.append(clean_name)

                    # An artificial residue will be a residue in a realchain but not included in force field map
                    if clean_name not in Resi_map2 and clean_name != "NAN":
                        self.if_art_resi = 1
                        ## PLACE HOLDER for further operation on artificial residue ##

            else:
                # Judge if containing any ligand
                for name in self.raw_sequence[chain]:
                    clean_name = name.strip(" ")
                    if clean_name not in TIP3P_map and clean_name != "NAN":
                        self.if_ligand = 1
                        ## PLACE HOLDER for further operation on artificial residue ##

            # only add realchain to the self.sequence
            if len(chain_seq) != 0:
                chain_Index = chr(new_index + 65)
                self.sequence[chain_Index] = chain_seq
                new_index = new_index + 1

    def _get_Oneletter(self):
        """
        (Used internally) convert sequences in self.sequence to oneletter-based str
        - The 'NAN' is convert to '-'
        - Present unnature residue as full 3-letter name

        save to self.sequence_one
        """
        if len(self.sequence) == 0:
            print("The self.sequence should be obtained first")
            raise IndexError

        for chain in self.sequence:
            chain_Seq = ""
            for name in self.sequence[chain]:
                c_name = name.strip(" ")
                if c_name == "NAN":
                    chain_Seq = chain_Seq + "-"
                else:
                    if c_name in Resi_map2:
                        chain_Seq = chain_Seq + Resi_map2[c_name]
                    else:
                        chain_Seq = chain_Seq + " " + c_name + " "

            self.sequence_one[chain] = chain_Seq

    def get_if_complete(self):
        """
        Judge if the self.sequence (from the get_seq) has internal missing parts
        Save the result to:
            self.if_complete ------- for intire PDB
            self.if_complete_chain - for each chain
        """
        if len(self.sequence.keys()) == 0:
            print("Please get the sequence first")
            raise IndexError

        self.if_complete = 1  # if not flow in then 1
        self.if_complete_chain = {}

        for chain in self.sequence:
            self.if_complete_chain[chain] = 1  # if not flow in then 1
            for resi in self.sequence[chain]:
                if resi == "NAN":
                    self.if_complete = 0
                    self.if_complete_chain[chain] = 0
                    break
        return self.if_complete, self.if_complete_chain

    def get_missing(self, seq):
        """
        get_missing(self, seq)
        Compare self.sequence (from the get_seq) with the seq (str from the uniport)
        1. No missing
        2. Terminal missing
        3. Internal missing
        """
        pass

    def PDB_loopmodel_refine(self, method="Rosetta"):
        """
        Use different methods to model the missing sequence
        Methods:
        - pyRosetta
        - trRosetta
        remember to update the name.
        """
        pass

    def get_protonation(
        self,
        ph: float = 7.0,
        ffout: str = "AMBER",
        keep_id: int = 0,
        if_prt_ligand: int = 1,
    ) -> None:
        # TODO check that the ph is in the range 0-14
        # self._get_protonation_pdb2pqr(ph=ph)
        # self._protonation_Fix(out_path, ph=ph, keep_id=keep_id, if_prt_ligand=if_prt_ligand)
        # self.path = out_path
        # self._update_name()
        # self.stru.name = self.name

        # def _get_protonation_pdb2pqr(self, ffout: str="AMBER", ph: float=7.0, out_path: str=""):
        """
        Use PDB2PQR to get the protonation state for current PDB. (self.path)
        current implementation just use the outer layer of PDB2PQR. Update to inner one and get more infomation in the furture. 
            (TARGET: 1. what is deleted from the structure // metal, ligand)
        
        save the result to self.pqr_path
        """
        # input of PDB2PQR
        pdb2pqr_parser = build_pdb2pqr_parser()
        args = pdb2pqr_parser.parse_args(
            [
                "--ff=PARSE",
                "--ffout=" + ffout,
                "--with-ph=" + str(ph),
                "--log-level=CRITICAL",
                self.path_name,
                self.pqr_path,
            ]
        )
        pdb2pqr_log = f"{self.work_dir}/pdb2pqr_output.log"
        _LOGGER.info(f"Running pdb2pqr on '{self.path_name}'...")
        run_pdb2pqr(args)
        _LOGGER.info(f"Finished running pdb2pqr! Output saved to '{self.pqr_path}'")
        # Add missing atom (from the PDB2PQR step. Update to func result after update the _get_protonation_pdb2pqr func)
        # Now metal and ligand
        old_stru = structure_from_pdb(self.no_water_path)
        new_stru = structure_from_pdb(self.pqr_path)
        # find Metal center and combine with the pqr file
        metal_list = old_stru.get_metal_center()
        if len(metal_list):
            _LOGGER.info(f"Merging {len(metal_list)} metal centers in old structure!")
            _LOGGER.info(f"Adding metal centers to new structure...")
            new_stru.add(metal_list, sort=0)
            _LOGGER.info(f"Metal centers added!")
            # fix metal environment
            _LOGGER.info(f"Protonating newly added metals...")
            new_stru.protonation_metal_fi(Fix=1)
            _LOGGER.info(f"Protonation complete!")

        # protonate ligands and combine with the pqr file
        ligand_list = old_stru.ligands
        if len(old_stru.ligands):
            _LOGGER.info(f"Merging {len(ligand_list)} ligands in old structure!")
            lig_dir = self.work_dir + "/ligands/"
            safe_mkdir(lig_dir)
            # TODO: logging
            new_ligands = list(
                map(
                    lambda ll: protonate_ligand(ll, dirname=lig_dir, ph=ph), ligand_list
                )
            )
            new_stru.add(new_ligands, sort=0)

        # PLACE HOLDER for other fix

        # build file
        if not keep_id:
            new_stru.sort()
        # TODO make this better for nameing
        new_stru.build(f"{self.work_dir}/final.pdb", keep_id=keep_id)
        self.stru = new_stru

    """
    ========
    Docking
    ========
    '''
    def Dock_Reactive_Substrate(self, substrate, reactive_define, local_lig = 0):
        '''
        Dock substrates into the apo enzyme in a reactive conformation
        Update self.path after docking. Will not update self.path in the middle of the docking
        -----
        In:     "apo-enzyme structure" self.path (pdb)
                "substrate structure" substrate (smile/sdf)
                "constraint parameters" reactive_define - a set of constraint information that defines *Reactive*
                    [reacting residue]
                        id
                        atom_name_1
                    [substrate]
                        atom_name_1
                    --- Will Generate ---
                    [reacting residue]
                        name 
                        atom_name_2
                        atom_name_3
                    [substrate]
                        id 
                        name 
                        atom_name_2 
                        atom_name_3
                    [distance]
                        ideal
                        allowed deviation
                        force constant
                        if_covalent
                    [angle]
                        ideal
                        allowed deviation
                        force constant
        Out:    Enzyme-Substrate Reactive Complex
        -----
        '''
        apo_enzyme = self.path
        if local_lig:
            lig_dir = self.dir+'/ligands/'
            met_dir = self.dir+'/metalcenters/'
        else:
            lig_dir = self.dir+'/../ligands/'
            met_dir = self.dir+'/../metalcenters/'
        mkdir(lig_dir)
        mkdir(met_dir)

        # Build ligands 
        cofactors = self.stru.build_ligands(lig_dir, ifcharge=1, ifunique=1)
        # cofactor_params, cofactors_pdb = generate_Rosetta_params(cofactors, lig_dir, resn='same', out_pdb_name='same')
        # # clean enzyme
        # self.fix_residue_names_for_rosetta()
        # # prepare substrate
        # sub_confs = Conformer_Gen(substrate, method='rdkit')
        # sub_params, sub_pdb, sub_confs_pdb = generate_Rosetta_params(sub_confs, lig_dir, resn='SUB', out_pdb_name='SUB')
        # # define reactive
        # -- search for RULE OF INPUT ID and get a map --
        # self.RDock_generate_cst_file(reactive_define, name_map?)        
        # # initial placement
        # self._RDock_initial_placement(sub_pdb)
        # self._RDock_add_remark_line()
        # # make config file for Rosetta
        # score_path = XXX
        # out_pdb_lib_path = XXX # control the output
        # option_path = Config.Rosetta.generate_xmlNoption(self.path, [cofactor_params, sub_params], score_path, nstruct=100, task='RDock')
        # # Run rosetta script
        # PDB.Run_Rosetta_Script(option_path)
        # Final_pdbs = self.RDock_get_result(score_path, out_pdb_lib_path)
        # self.Rosetta_pdb_to_Amber()



    '''
    ========
    Mutation
    ========
    """

    # def PDB2PDBwLeap(self):
    def apply_mutations(self) -> None:
        """
        Apply mutations using tleap. Save mutated structure PDB in self.path
        ------------------------------
        Use MutaFlag in self.MutaFlags
        Grammer (from Add_MutaFlag):
        X : Original residue name. Leave X if unknow. 
            Only used for build filenames. **Do not affect any calculation.**
        A : Chain index. Determine by 'TER' marks in the PDB file. (Do not consider chain_indexs in the original file.)
        11: Residue index. Strictly correponding residue indexes in the original file. (NO sort applied)
        Y : Target residue name.  

        **WARNING** if there are multiple mutations on the same index, only the first one will be used.
        """

        # Prepare a label for the filename
        tot_Flag_name = "_".join(list(map(mutaflag_to_str, self.mutations)))
        # Operate the PDB
        out_PDB_path1 = (
            self.work_dir + "/" + self.base_pdb_name + tot_Flag_name + "_tmp.pdb"
        )
        out_PDB_path2 = self.path_name + tot_Flag_name + ".pdb"
        chain_count = 1
        pdb_lines = read_pdb_lines(self.path_name)
        mask = [True] * len(pdb_lines)
        for pdb_l in pdb_lines:
            if pdb_l.is_TER():
                chain_count += 1
            match = 0
            # only match in the dataline and keep all non data lines
            if not pdb_l.is_ATOM():
                continue
            for mf in self.mutations:
                # Test for every Flag for every lines

                if (
                    chr(64 + chain_count) == mf.chain_index
                    and pdb_l.resi_id == mf.residue_index
                ):
                    # do not write old line if match a MutaFlag
                    match = 1
                    # Keep OldAtoms of targeted old residue
                    target_residue = mf.target_residue
                    # fix for mutations of Gly & Pro
                    old_atoms = {
                        "G": ["N", "H", "CA", "C", "O"],
                        "P": ["N", "CA", "HA", "CB", "C", "O"],
                    }.get(mf.target_residue, ["N", "H", "CA", "HA", "CB", "C", "O"])

                    line = pdb_l.line
                    for oa in old_atoms:
                        if oa == pdb_l.atom_name:
                            pdb_l.line = f"{line[:17]}{convert_to_three_letter(mf.target_residue)}{line[20:]}"

        write_lines(out_PDB_path1, list(map(lambda pl: pl.line, pdb_lines)))
        leapin_path = f"{self.work_dir}/leap_P2PwL.in"
        leap_lines = [
            "source leaprc.protein.ff14SB",
            f"a = loadpdb {out_PDB_path1}",
            f"savepdb a {out_PDB_path2}",
            "quit",
        ]
        write_lines(leapin_path, leap_lines)
        em.run_command(
            "tleap", [f"-s -f {leapin_path} > {self.work_dir}/leap_P2PwL.out"]
        )
        safe_rm("leap.log")

    def generate_mutations(
        self,
        n: int,
        muta_flags: List[MutaFlag] = None,
        restrictions: str = None,
        random_state=100,
    ) -> None:
        # TODO add restrictions
        # TODO check for restrictions
        # TODO setup to be compliant with the existing API grammar
        existing = set()
        temp = []
        if muta_flags:
            # TODO convert and check if the mutations are valid
            # TODO check the target residue exists in the structure
            for mf in muta_flags:
                key = (mf.chain_index, mf.residue_index)
                if key in existing:
                    _LOGGER.warn(f"Multiple mutations supplied at location {key}")
                else:
                    existing.add(key)
                    temp.append(mf)

        muta_flags = temp

        if muta_flags and len(muta_flags) >= n:
            _LOGGER.warn(
                f"Supplied mutation flags meet or exceed the number of desired mutations"
            )
            self.mutations = muta_flags
            return

        candidates = self.stru.all_possible_mutations()
        np.random.seed(random_state)
        np.random.shuffle(candidates)

        while len(muta_flags) < n and len(candidates):
            curr = candidates.pop()
            key = (curr.chain_index, curr.residue_index)
            if key in existing:
                continue
            else:
                muta_flags.append(curr)
                existing.add(key)

        if len(muta_flags) != n:
            _LOGGER.warn(
                f"Unable to generate enough mutations. Missing {n-len(muta_flags)}"
            )

        self.mutations = muta_flags
        assert len(self.mutations) == len(set(self.mutations))

    def Add_MutaFlag(self, Flag: str = "r", if_U: bool = 0, if_self: bool = 0):
        """Determine which mutation to deploy to the structure.
        
        User can 1. assign specific mutation(s). 
              or 2. assign random mutation(s) by specifing a rule of randomlization.
        The function will 
        1. decode the input to a list of python tuples that 
            each contains the information of a mutation as ('X','A','###','Y') 
            save to *self.MutaFlags*
        2. check the if the assigned mutation is reasonable that
            specifing the correct original residue
            within the range of chain ids
            within the range of residue ids
            the target residue is within the range of canonical AAs

        Args:
        (self):
        The requirements for the input pdb object are:
            TODO(shaoqz): add after refactoring
        Flag: 
            str or list of strings indicating specifc mutation(s) or a rule of randomlization. (default: r)
            Grammer:
            **Assign specific mutation(s)**
                'XA##Y'
                For each str that represent a mutation, the format should be 'XA##Y', in which
                X : Original residue letter (One-letter code). Leave X if unknow. 
                    (Only used for checking that will pop a warning if the residue doesn't match, 
                    which usually means you are using a wrong residue id and requires a double check.)
                A : Chain index. This is optional witha default value: A.
                    (Different chains are defined by 'TER' marks in the PDB file and the index is noted in order.
                    * Note this will not always be the same as the chain id in the PDB line
                    * defined in self.stru)
                ##: Residue index. 
                    (Strictly correponding residue indexes in the original file.
                    * defined in self.stru)
                Y : Target residue letter.
                    (The residue type after the mutation.) 
            example:
                >>> pdb_obj.Add_MutaFlag('D83K')
                >>> pdb_obj.MutaFlags
                [('D', 'A', '83', 'K')]
                >>> pdb_obj.Add_MutaFlag('DA83K')
                >>> pdb_obj.MutaFlags
                [('D', 'A', '83', 'K')]
                >>> pdb_obj.Add_MutaFlag(['DA83K', 'EB226P'])
                >>> pdb_obj.MutaFlags
                [('D', 'A', '83', 'K'), ('E', 'B', '226', 'P')]
            
            **Assign random mutation(s)** 
                ['r', {'position': 'keyword', 
                       'target_resi':'keyword'}]
                To assigne random mutaion, start the list with 'r' or 'random' followed by a map that
                defines the rule of the randomlization. In the map dictionary there are 2 keys to fill:
                position        : use a keyword or a pattern to define availiable positions to mutate.
                target_residue  : use a keyword or a pattern to define availiable target residues.
                (* Note that when using random assignment, the list can no longer contain any manual assign str for mutation.
                e.g.: pdb_obj.Add_MutaFlag(['V23T', 'r']) is not valid)
            example:
                TODO(shaoqz): finish this function.
        if_U: 
            if include mutations to U (selenocysteine) in random generation.
        if_self:
            if "mutation to the same amino acid" is allowed for random mutation.

        Returns:
        There are two parts of actually outputs:
        In self.MutaFlags: Assigned mutations in a list of tuples
            example: 
            [('D', 'A', '83', 'K'), ('E', 'B', '226', 'P')]
        In return value: A str representing a tag of current mutations to the PDB.
            example:
            _DA83K_EB226P
        """

        if type(Flag) == str:

            if Flag == "r" or Flag == "random":
                resi_1 = ""
                resi_2 = ""
                Muta_c_id = ""
                Muta_r_id = ""

                # Flag Generation：
                # Random over the self.stru. Strictly correponding residue indexes in the original file. (no sort applied)
                self.get_stru()
                # random over the structure "protein" part.
                chain = choice(self.stru.chains)
                resi = choice(chain.residues)
                Muta_c_id = chain.id
                Muta_r_id = str(resi.id)
                if resi.name in Resi_map2:
                    resi_1 = Resi_map2[resi.name]
                else:
                    if Config.debug >= 1:
                        print(
                            "WARNING: pdb.Add_MutaFlag(): A non-canonical animo acid is being mutated! The MutaFlag will have a 3-letter code for the original residue ."
                        )
                    resi_1 = resi.name
                # random over the residue list
                if if_U:
                    m_Resi_list = Resi_list
                else:
                    m_Resi_list = Resi_list[:-1]
                resi_2 = choice(m_Resi_list)
                # avoid or not mutation to self
                if not if_self:
                    while resi_2 == resi_1:
                        resi_2 = choice(m_Resi_list)

                MutaFlag = (resi_1, Muta_c_id, Muta_r_id, resi_2)
                self.MutaFlags.append(MutaFlag)

            else:
                MutaFlag = self._read_MutaFlag(Flag)
                self.MutaFlags.append(MutaFlag)

        if type(Flag) == list:
            for i in Flag:
                MutaFlag = self._read_MutaFlag(i)
                self.MutaFlags.append(MutaFlag)

        if Config.debug >= 1:
            print("Current MutaFlags:")
            for flag in self.MutaFlags:
                print(self._build_MutaName(flag))

        label = ""
        for flag in self.MutaFlags:
            label = label + "_" + self._build_MutaName(flag)
        return label

    def _read_MutaFlag(self, Flag):
        """
        decode the manually input MutaFlag. Return (resi_1, chain_id, resi_id, resi_2)
        --------------
        Grammer(XA11Y):
        X : Original residue name. Leave X if unknow. 
            Only used for build filenames. **Do not affect any calculation.**
        A : Chain index. Determine by 'TER' marks in the PDB file. (Do not consider chain_indexs in the original file.)
        11: Residue index. Strictly correponding residue indexes in the original file. (NO sort applied)
        Y : Target residue name.

        the later 3 will go through a san check to make sure they are within the range:
        A : chain.id range in self.stru.chains
        11: resi.id range in self.stru.chain[int].residues
        Y : Resi_list
        """
        pattern = r"([A-Z])([A-Z])?([0-9]+)([A-Z])"
        F_match = re.match(pattern, Flag)
        if F_match is None:
            raise Exception(
                "_read_MutaFlag: Required format: XA123Y (or X123Y indicating the first chain)"
            )

        resi_1 = F_match.group(1)
        chain_id = F_match.group(2)
        resi_id = str(F_match.group(3))
        resi_2 = F_match.group(4)

        # default
        if F_match.group(2) is None:
            chain_id = "A"
            if Config.debug >= 1:
                print(
                    "_read_MutaFlag: No chain_id is provided! Mutate in the first chain by default. Input: "
                    + Flag
                )

        # san check of the manual input
        self.get_stru()
        chain_id_list = [i.id for i in self.stru.chains]
        if not chain_id in chain_id_list:
            raise Exception(
                "_read_MutaFlag: San check failed. Input chain id in not in range."
                + line_feed
                + " range: "
                + repr(chain_id_list)
            )
        chain_int = ord(chain_id) - 65
        resi_id_list = [str(i.id) for i in self.stru.chains[chain_int].residues]
        if not resi_id in resi_id_list:
            raise Exception(
                "_read_MutaFlag: San check failed. Input resi id in not in range."
                + line_feed
                + " range: "
                + repr(resi_id_list)
            )
        if not resi_2 in Resi_list:
            raise Exception(
                "_read_MutaFlag: Only support mutate to the known 21 residues. AmberMaps.Resi_list: "
                + repr(Resi_list)
            )

        return (resi_1, chain_id, resi_id, resi_2)

    def _build_MutaName(self, Flag):
        """
        Take a MutaFlag Tuple and return a str of name
        """
        return Flag[0] + Flag[1] + Flag[2] + Flag[3]

    def rm_allH(self, ff="Amber", ligand=False):
        # TODO: CJ: should this get called by mutation method?
        """
        remove wrong hydrogens added by leap after mutation. (In the case that the input file was a H-less one from crystal.)
        ----------
        if_ligand:
        0 - remove Hs of standard protein residues only.
        1 - remove all Hs base on the nomenclature. (start with H and not in the non_H_list)
        """
        # out path
        o_path = self.path_name + "_rmH.pdb"
        pdb_lines: List[PDBLine] = read_pdb_lines(self.curr_path)
        mask = [True] * len(pdb_lines)
        # crude judgement of H including customized H
        if ligand:
            not_H_list = ["HG", "HF", "HS"]  # non-H elements that start with "H"
            for idx, pl in enumerate(pdb_lines):
                if pl.is_ATOM():
                    atom_name = line[12:16].strip()
                    if (
                        pl.atom_name.startswith("H")
                        and pl.atom_name[:2] not in not_H_list
                    ):
                        mask[idx] = False
        else:
            H_aliases = get_element_aliases(ff, "H")
            for idx, pl in enumerate(pdb_lines):
                if pl.atom_name in H_aliases and pl.is_residue_line():
                    mask[idx] = False

        pdb_lines = np.array(pdb_lines)[mask]
        write_lines(o_path, list(map(str, pdb_lines)))

    """
    ========
    General MD
    ========
    """

    @classmethod
    def get_charge_list(cls, prmtop_path):
        """
        Get charge from the .prmtop file
        Take the (path) of .prmtop file and return a [list of charges] with corresponding to the atom sequence
        -----------------
        * Unit transfer in prmtop: http://ambermd.org/Questions/units.html
        """
        with open(prmtop_path) as f:

            charge_list = []
            line_index = 0

            for line in f:

                line_index = line_index + 1  # current line

                if line.strip() == r"%FLAG POINTERS":
                    format_flag = line_index
                if line.strip() == r"%FLAG CHARGE":
                    charge_flag = line_index

                if "format_flag" in dir():
                    if line_index == format_flag + 2:
                        N_atom = int(line.split()[0])
                        del format_flag

                if "charge_flag" in dir():
                    if (
                        line_index >= charge_flag + 2
                        and line_index <= charge_flag + 1 + ceil(N_atom / 5)
                    ):
                        for i in line.strip().split():
                            charge_list.append(float(i) / 18.2223)
        return charge_list

    """
    ========
    QM Cluster
    ========    
    """

    def PDB2QMCluster(
        self,
        atom_mask,
        spin=1,
        o_dir="",
        tag="",
        QM="g16",
        g_route=None,
        ifchk=0,
        val_fix="internal",
    ):
        """
        Build & Run QM cluster input from self.mdcrd with selected atoms according to atom_mask
        ---------
        spin: specific spin state for the qm cluster. (default: 1)
        QM: QM engine (default: Gaussian)
            g_route: Gaussian route line
            ifchk: if save chk file of a gaussian job. (default: 0)
        val_fix: fix free valance of truncated strutures
            = internal: add H to where the original connecting atom is. 
                        special fix for classical case:
                        "sele by residue" (cut N-C) -- adjust dihedral for added H on N.
            = openbabel TODO
        ---data---
        self.frames
        self.qm_cluster_map (PDB atom id -> QM atom id)
        """
        # make folder
        if o_dir == "":
            o_dir = self.dir + "/QM_cluster" + tag
        mkdir(o_dir)
        # update stru
        self.get_stru()
        # get sele
        if val_fix == "internal":
            sele_lines, sele_map = self.stru.get_sele_list(
                atom_mask, fix_end="H", prepi_path=self.prepi_path
            )
        else:
            sele_lines, sele_map = self.stru.get_sele_list(atom_mask, fix_end=None)
        self.qm_cluster_map = sele_map
        # get chrgspin
        chrgspin = self._get_qmcluster_chrgspin(sele_lines, spin=spin)
        if Config.debug >= 1:
            print("Charge: " + str(chrgspin[0]) + " Spin: " + str(chrgspin[1]))

        # make inp files
        frames = Frame.fromMDCrd(self.mdcrd)
        self.frames = frames
        if QM in ["g16", "g09"]:
            gjf_paths = []
            if Config.debug >= 1:
                print("Writing QMcluster gjfs.")
            for i, frame in enumerate(frames):
                gjf_path = o_dir + "/qm_cluster_" + str(i) + ".gjf"
                frame.write_sele_lines(
                    sele_lines,
                    out_path=gjf_path,
                    g_route=g_route,
                    chrgspin=chrgspin,
                    ifchk=ifchk,
                )
                gjf_paths.append(gjf_path)
            # Run inp files
            qm_cluster_out_paths = PDB.Run_QM(gjf_paths, prog=QM)
            # get chk files if ifchk
            if ifchk:
                qm_cluster_chk_paths = []
                for gjf in gjf_paths:
                    chk_path = gjf[:-3] + "chk"
                    qm_cluster_chk_paths.append(chk_path)
        if QM == "ORCA":
            pass

        self.qm_cluster_out = qm_cluster_out_paths
        if ifchk:
            self.qm_cluster_chk = qm_cluster_chk_paths
            return self.qm_cluster_out, self.qm_cluster_chk

        return self.qm_cluster_out

    def _get_qmcluster_chrgspin(self, sele, spin=1):
        """
        get charge for qmcluster of sele from self.prmtop
        """
        # get chrg list
        chrg_list_all = PDB.get_charge_list(self.prmtop_path)
        # sum with sele
        sele_chrg = 0
        for sele_atom in sele.keys():
            if "-" in sele_atom:
                continue
            # skip backbone atoms (cut CA-CB when calculate charge)
            if sele_atom[-1] == "b":
                continue

            # clean up
            if sele_atom[-1] not in "1234567890":
                sele_id = sele_atom[:-1]
            else:
                sele_id = sele_atom
            sele_chrg += chrg_list_all[int(sele_id) - 1]

        return (round(sele_chrg), spin)

    @classmethod
    def Run_QM(cls, inp, prog="g16"):
        """
        Run QM with prog
        """
        if prog == "g16":
            outs = []
            for gjf in inp:
                out = gjf[:-3] + "out"
                if Config.debug > 1:
                    print(
                        "running: "
                        + Config.Gaussian.g16_exe
                        + " < "
                        + gjf
                        + " > "
                        + out
                    )
                os.system(Config.Gaussian.g16_exe + " < " + gjf + " > " + out)
                outs.append(out)
            return outs

        if prog == "g09":
            outs = []
            for gjf in inp:
                out = gjf[:-3] + "out"
                if Config.debug > 1:
                    print(
                        "running: "
                        + Config.Gaussian.g09_exe
                        + " < "
                        + gjf
                        + " > "
                        + out
                    )
                os.system(Config.Gaussian.g09_exe + " < " + gjf + " > " + out)
                outs.append(out)
            return outs

    """
    ========
    QM Analysis 
    ========
    """

    def get_field_strength(
        self, atom_mask, a1=None, a2=None, bond_p1="center", p1=None, p2=None, d1=None
    ):
        """
        use frame coordinate from *mdcrd* and MM charge from *prmtop* to calculate the field strength of *p1* along *p2-p1* or *d1*
        atoms in *atom_mask* is included. (TODO: or an exclude one?)
        -------------------------------------
        a1 a2:  id of atoms compose the bond
        bond_p1:method to generate p1
                - center
                - a1
                - TODO exact dipole center of the current structure
                - ...
        p1:     the point where E is calculated
        p2:     a point to fix d1
        d1:     the direction E is projected
        return an ensemble field strengths
        """
        Es = []
        chrg_list = PDB.get_charge_list(self.prmtop_path)

        # san check
        if a1 == None and p1 == None:
            raise Exception(
                "Please provide a 1nd atom (a1=...) or point (p1=...) where E is calculated "
            )
        if a2 == None and p2 == None and d1 == None:
            raise Exception(
                "Please provide a 2nd atom (a2=...) or point (p2=...) or a direction (d1=...)"
            )
        if a1 == None and a2 == None and bond_p1 == "center":
            raise Exception(
                "Please provide a both atom (a1=... a2=...) when bond_p1 = center; Or change bond_p1 to a1 to calculate E at a1"
            )
        if a1 != None and a2 != None and bond_p1 not in ["center", "a1"]:
            raise Exception("Only support p1 selection in center or a1 now")

        if self.frames == None:
            self.frames = Frame.fromMDCrd(self.mdcrd)

        # decode atom mask (stru corresponding to mdcrd structures)
        atom_list = decode_atom_mask(self.stru, atom_mask)

        for frame in self.frames:
            # get p2
            if a2 != None:
                p2 = frame.coord[a2 - 1]
            # get p1
            if a1 != None:
                if bond_p1 == "a1":
                    p1 = frame.coord[a1 - 1]
                if bond_p1 == "center":
                    p1 = get_center(frame.coord[a1 - 1], p2)
                if bond_p1 == "xxx":
                    pass
            # sum up field strength
            E = 0
            for atom_id in atom_list:
                # search for coord and chrg
                coord = frame.coord[atom_id - 1]
                chrg = chrg_list[atom_id - 1]
                E += get_field_strength_value(coord, chrg, p1, p2=p2, d1=d1)
            Es.append(E)

        return Es


def get_PDB(name):
    """
    connect to the database
    """
    pass


def PDB_to_AMBER_PDB(path):
    """
    Make the file convertable with tleap without error
    - Test result: The header part cause duplication of the residues. Deleting that part may give normal tleap output
    - Test result: For some reason, some ligand will miss if directly covert after simple header cutting
    - Test result: `cat 1NVG.pdb |grep "^ATOM\\|^HETATM\\|^TER|^END" > 1NVG-grep.pdb` will be fine
    - WARNINGs
    """
    pass