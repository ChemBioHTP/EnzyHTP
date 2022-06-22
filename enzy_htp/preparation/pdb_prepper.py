"""Defines a PDBPrepper class that takes in a raw pdb and prepares it for analysis.
For existing users of EnzyHTP, this largely replaces what Class_PDB did.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-12
"""
import shutil
import numpy as np

import pdb2pqr
from typing import Set, Union, List, Tuple
from enzy_htp.core import _LOGGER

from .protonate import protonate_pdb, protonate_missing_elements, check_valid_ph

import enzy_htp.core as core

from enzy_htp.core import file_system as fs
from enzy_htp.structure import (
    compare_structures,
    merge_right,
    structure_from_pdb,
    Ligand,
    ligand_from_pdb,
    Chain,
    Structure,
)

from .pdb_line import PDBLine, read_pdb_lines


#from .mutate import (
#    Mutation,
#    mutation_to_str,
#    decode_mutaflags,
#    get_all_combinations,
#    get_all_combinations,
#)
#
from enzy_htp.chemical import convert_to_three_letter, get_element_aliases


class PDBPrepper:
    """Class that handles initial preparation of a specified Structure().

    Attributes:
    current_path_ : The path to the most recently saved version of the Structure().
        work_dir : The scratch directory where temporary files are saved. The current time in the format "YYYY_MM_DD" is used if the keyword is not specified.
            all_paths : A list of all paths created by the PDBPrepper() object. Stored in reverse order of age.

    """

    def __init__(self, pdb_name, **kwargs):
        """Inits PDBPrepper with a pdb file and optionally a work directory to place temporary files."""
        self.current_path_: str = None
        self.no_water_path = None
        self.pdb_path = pdb_name
        self.base_pdb_name = fs.base_file_name(pdb_name)
        self.work_dir = kwargs.get(
            "work_dir", f"{fs.get_current_time()}_{self.base_pdb_name}"
        )
        fs.safe_mkdir(self.work_dir)
        shutil.copy(self.pdb_path, f"{self.work_dir}/{self.base_pdb_name}.pdb")
        self.path_name = f"{self.work_dir}/{self.base_pdb_name}.pdb"
        self.pqr_path = str()
        #self.mutations = []
        self.current_path_ = self.path_name
        self.all_paths = [self.current_path_]

    def current_path(self) -> str:
        """Gets the most recent PDB that has been created."""
        return self.current_path_

    def rm_wat(self) -> str:
        """Removes water and ions from the PDB and saves a version in the specified work_dir. Returns the path to the non water PDB."""
        pdb_lines: List[PDBLine] = read_pdb_lines(self.path_name)
        pdb_lines = list(
            filter(lambda pl: not (pl.is_water() or pl.is_CRYST1()), pdb_lines)
        )
        mask = [True] * len(pdb_lines)

        for pidx, pl in enumerate(pdb_lines[:-1]):
            if pl.is_TER() and pdb_lines[pidx + 1].is_TER():
                mask[pidx + 1] = False

        pdb_lines = np.array(pdb_lines)[mask]
        self.no_water_path = f"{self.work_dir}/{self.base_pdb_name}_rmW.pdb"
        fs.write_lines(self.no_water_path, list(map(str, pdb_lines)))

        self.current_path_ = self.no_water_path
        self.all_paths.append(self.current_path_)
        return self.current_path_

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
                            Chain_sequence = (
                                Chain_sequence
                                + [
                                    "NAN",
                                ]
                                * missing_length
                            )

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

        pass

    def get_protonation(
        self,
        ph: float = 7.0,
        ffout: str = "AMBER",
    ) -> str:
        """Uses PDB2PQR to get the protonation state for a PDB file, saving the protonated .pdb file.
        Separately runs PDB2PQR on the ligand and metal and combines the output. Returns path to the protonated
        PDB file."""
        # input of PDB2PQR
        check_valid_ph(ph)
        self.pqr_path = f"{fs.remove_ext(self.current_path_)}.pqr.pdb"
        self.all_paths.append(self.pqr_path)
        protonate_pdb(self.path_name, self.pqr_path)
        # Add missing atom (from the PDB2PQR step. Update to func result after update the _get_protonation_pdb2pqr func)
        # Now metal and ligand
        new_structure: Structure = protonate_missing_elements(
            self.no_water_path, self.pqr_path, self.work_dir
        )

        self.current_path_ = f"{fs.remove_ext(fs.remove_ext(self.pqr_path))}_aH.pdb"
        self.all_paths.append(self.current_path_)

        new_structure.to_pdb(self.current_path_)
        self.stru = new_structure

        return self.current_path_

    """
    ========
    Mutation
    ========
    """



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
