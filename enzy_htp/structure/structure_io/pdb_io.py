"""Generation/construction of Structure objects from PDB files and exporting it vice versa
Definition of PDB file format (http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)
Also contains util function that requires knoweledge of PDB format

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2022-08-01
"""
from collections import defaultdict
import os
import string
import sys
from typing import Dict, List, Tuple, Union
from plum import dispatch
from biopandas.pdb import PandasPdb
import pandas as pd

from enzy_htp.core import _LOGGER
from enzy_htp.core.exception import IndexMappingError
import enzy_htp.core.file_system as fs
from enzy_htp.core.general import if_list_contain_repeating_element, list_remove_adjacent_duplicates
from enzy_htp.core.pandas_helper import batch_edit_df_loc_value, split_df_base_on_column_value
import enzy_htp.chemical as chem
from ._interface import StructureParserInterface
from ..atom import Atom
from ..metal_atom import MetalUnit, residue_to_metal
from ..residue import Residue
from ..solvent import Solvent, residue_to_solvent
from ..ligand import Ligand, residue_to_ligand
from ..chain import Chain
from ..structure import Structure


class PDBParser(StructureParserInterface):
    """
    Parser for covert PDB into Structure and vice versa
    """

    def __init__(self) -> None:  # pylint: disable=super-init-not-called
        """
        not sure what can be init here. Maybe some configuration?
        """
        pass

    # knowledge
    PDB_NONSTRU_INFO_LINE_NAME = [
        "HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "NUMMDL", "AUTHOR", "REVDAT", "SPRSDE", "JRNL", "REMARK", "DBREF",
        "SEQADV", "SEQRES", "HELIX", "SHEET", "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3"
    ]  #these keywords are for general information in the PDB database

    # interface
    @classmethod
    def get_structure(cls,
                      path: str,
                      model: int = 0,
                      add_solvent_list: List = None,
                      add_ligand_list: List = None,
                      remove_trash: bool = True,
                      give_idx_map: bool = False,
                      allow_multichain_in_atom: bool = False) -> Union[Structure, tuple]:
        """
        Converting a PDB file (as its path) into the Structure()
        Arg:
            path:
                the file path of the PDB file
            model:
                The selected model index if there are multiple models in the file
                (start from 0 default: 0)
                (assume MODEL appears in order in the file)
            add_solvent_list:
                (used for categorize residues) additional names for solvent
                solvents are recognized by matching names in these lists in
                non-polypeptide chains:
                    chem.RD_SOLVENT_LIST + add_solvent_list
            add_ligand_list:
                (used for categorize residues) additional names for ligands
                ligands are recognized by not matching names in these lists
                in non-polypeptide chains:
                    1. chem.RD_SOLVENT_LIST + add_solvent_list
                    2. chem.RD_NON_LIGAND_LIST - add_ligand_list
                    * solvent list have higher pirority
            remove_trash:
                if remove trash ligands
            give_idx_map:
                if return a tuple of (Structure, idx_change_mapper) - which is a dict 
                with {(old_chain_id, old_residue_id): (new_chain_id, new_residue_id), ... }
            allow_multichain_in_atom:
                (used for resolving chain id)
                if allowing multiple chain IDs appears in the same chain that consists of ATOM
                records. It conflict with the definiation of PDB file format but it is useful
                for resolving chain id of pymol2 exported multichain PDBs.

        Return:
            Structure()
        """
        _LOGGER.debug(f"working on {path}")
        cls._check_valid_pdb(path)
        # covert to dataframe with biopanda (maybe write our own in the future)
        input_pdb = PandasPdb()
        input_pdb.read_pdb(path)

        #region (PDB conundrums)
        # deal with multiple model
        target_model_df, target_model_ter_df = cls._get_target_model(input_pdb.df, model=model)

        # TODO address when atom/residue/chain id is disordered respect to the line index
        # this is important to be aligned since Amber will force them to align and erase the chain id
        # a workaround is to add them back after Amber process
        idx_change_mapper = cls._resolve_missing_chain_id(
            target_model_df, target_model_ter_df, allow_multichain_in_atom=allow_multichain_in_atom)  # add missing chain id in place
        cls._resolve_alt_loc(target_model_df)  # resolve alt loc record in place by delele redundant df rows
        #endregion (Add here to address more problem related to the PDB file format here)

        # target_model_df below should be "standard" --> start building
        # mapper is for indicating superior information of the current level: e.g.: indicate
        # which residue and chain the atom belongs to in atom_mapper. Another workaround in
        # old enzyhtp is build iteratively chain frist and residues are build in the chain
        # builder and so on for atom. But current methods is used for better readability
        atom_mapper: Dict[tuple, Atom] = cls._build_atoms(target_model_df)
        res_mapper = cls._build_residues(atom_mapper, add_solvent_list, add_ligand_list)
        chain_list: Dict[str, Chain] = cls._build_chains(res_mapper, remove_trash)

        if give_idx_map:
            return Structure(chain_list), idx_change_mapper
        return Structure(chain_list)

    @classmethod
    @dispatch
    def get_file_str(cls, stru: Chain, if_renumber: bool = True, if_fix_atomname: bool = True) -> str:  # pylint: disable=function-redefined
        """
        dispatch for supporting get pdb file str with Chain only
        """
        if if_renumber:
            stru.parent.renumber_atoms()
        result_str = cls._write_pdb_chain(stru)
        result_str += f"END{os.linesep}"
        return result_str

    @classmethod
    @dispatch
    def get_file_str(cls, stru: Residue, if_renumber: bool = True, if_fix_atomname: bool = True) -> str:  # pylint: disable=function-redefined
        """
        dispatch for supporting get pdb file str with Residue only
        """
        if if_renumber:
            stru.parent.parent.renumber_atoms()
        result_str = cls._write_pdb_residue(stru)
        result_str += f"TER{os.linesep}END{os.linesep}"
        return result_str

    @classmethod
    @dispatch
    def get_file_str(cls, stru: Atom) -> str:  # pylint: disable=function-redefined
        """
        dispatch for supporting get pdb file str with Atom only
        """
        result_str = cls._write_pdb_atom(stru)
        result_str += f"TER{os.linesep}END{os.linesep}"
        return result_str

    @classmethod
    @dispatch
    def get_file_str(
            cls,
            stru: Structure,  # pylint: disable=function-redefined
            if_renumber: bool = True,
            if_fix_atomname: bool = True) -> str:
        """
        Convert Structure() into PDB file string. Only the simplest function is need for
        enzyme modeling.
        TODO support fixing all atom names before writing
        TODO do we need to add a mode where all ligand,metal,solvent are written into seperate 
        chains? we if encounter any need
            A better way to do this is to make derivate Parser classes and overwrite this as needed
            or use things like pdb4amber to accommodate software as needed.
        Return:
            a string of the PDB file
        """
        stru.sort_chains()
        if if_renumber:
            stru.renumber_atoms()
        result_str = ""
        for chain in stru:
            chain: Chain
            result_str += cls._write_pdb_chain(chain)
        result_str += f"END{os.linesep}"
        return result_str

    #region == pdb -> Stru ==
    @staticmethod
    def _check_valid_pdb(pdb_path: str) -> None:
        """Helper function that ensures the supplied pdb_path contain a valid pdb file.
        Private to structure.structure_io.pdb_io.py. Should NOT be called externally."""
        # check for right extension
        ext: str = fs.get_file_ext(pdb_path)
        if ext.lower() != ".pdb":
            _LOGGER.error(f"Supplied file '{pdb_path}' is NOT a PDB file. Exiting...")
            sys.exit(1)
        #check for file existance
        if not os.path.exists(pdb_path):
            _LOGGER.error(f"Supplied file '{pdb_path}' does NOT exist. Exiting...")
            sys.exit(1)
        # check for character encoding
        for idx, ll in enumerate(fs.lines_from_file(pdb_path)):
            if not ll.isascii():
                _LOGGER.error(f"The PDB '{pdb_path}' contains non-ASCII text and is invalid in line {idx}: '{ll}'. Exiting...")
                sys.exit(1)

    @staticmethod
    def _get_target_model(df: pd.DataFrame, model: int) -> Union[pd.DataFrame, None]:
        """
        do nothing when there is only one model in the pdb file (not MODEL line)
        get the #model model part of the dataframe when there are multiple models
        in the pdb file.
        Arg:
            df: the dataframe of a pdb file
            model: the index of target model
                (start from 0)
                (assume MODEL appears in order in the file)
        Return: (target_mdl_df, target_mdl_ter_df)
            target_mdl_df: return a copy of Atom and HETATM records from the input
                           dataframe correponding to the target model
            target_mdl_ter_df: return a copy of TER records from the input
                           dataframe correponding to the target model
        """
        if "MODEL" in df["OTHERS"]["record_name"].values:
            # san check for start/end indicator
            mdl_ind_sorted_w_lines = df["OTHERS"][(df["OTHERS"].record_name == "MODEL") |
                                                  (df["OTHERS"].record_name == "ENDMDL")].sort_values(by="line_idx", inplace=False)
            # check if they appear alternately
            first = 1
            for ind in mdl_ind_sorted_w_lines.values:
                if first:
                    last_ind = ind
                    first = 0
                    continue
                if ind[0] == last_ind[0]:
                    _LOGGER.error(f"Unclosed MODEL section is detected in pdb file. line: {last_ind[-1]}")
                    sys.exit(1)
                last_ind = ind
            # get target model unit range
            mdl_start_lines = list(df["OTHERS"][df["OTHERS"].record_name == "MODEL"]["line_idx"])
            mdl_end_lines = list(df["OTHERS"][df["OTHERS"].record_name == "ENDMDL"]["line_idx"])
            mdl_range = list(zip(mdl_start_lines, mdl_end_lines))
            target_mdl_range = mdl_range[model]
            mdl_query_pattern = f"line_idx > {target_mdl_range[0]} & line_idx < {target_mdl_range[1]}"
            # get model dataframe section as a copy
            target_mdl_df = pd.concat((df["ATOM"], df["HETATM"]), ignore_index=True).query(mdl_query_pattern).copy()
            target_mdl_ter_df = df["OTHERS"][df["OTHERS"].record_name == "TER"].query(mdl_query_pattern).copy()
        else:
            # get all dataframe as a copy if there"s no MODEL record
            target_mdl_df = pd.concat((df["ATOM"], df["HETATM"]), ignore_index=True).copy()  # insure unique df id
            target_mdl_ter_df = df["OTHERS"][df["OTHERS"].record_name == "TER"].copy()

        return target_mdl_df, target_mdl_ter_df

    @classmethod
    def _resolve_missing_chain_id(cls, df: pd.DataFrame, ter_df: pd.DataFrame, allow_multichain_in_atom: bool = False) -> None:
        """
        Function takes the dataframes of chains and ensures consistent naming
        of chains with no blanks.

        1. add missing chain id into the Atomic "df" when there is none
        2. update repeating chain id when it is seperated by "TER" 
        (e.g. chain A TER chain A -(update)-> chain A TER chain B. see 1Q4T.)
        according to the TER record provide by "ter_df"
        * ter df should share the same line index references as df (from same file)
        3. allow different ligand in the same chain
        ==details==
        There are in total 12 cases
            for each chain:
            contain atom missing chain id:
                have 0 chain id total:
                    add chain id
                have 1 chain id total:
                    repeat:
                        HET chain:
                            update to a new chain id for whole chain
                        ATOM chain:
                            abort
                    not repeat:
                        add the same chain id to atoms that missing
                have more than 1 chain id total:
                    abort
            contain no atom missing chain id:
                have 1 chain id total:
                    HET:
                        repeat:
                            update to a new chain id
                        not repeat:
                            record
                    ATOM:
                        repeat:
                            abort
                        not repeat:
                            record
                have more than 1 chain id total:
                    HET:
                        treat as individual chain. for each:
                        repeat:
                            update
                        not repeat:
                            record
                    ATOM:
                        not allow_multichain_in_atom (default)
                            abort
                        allow_multichain_in_atom
                            same as HET case. treat as individual chain.
                            repeat:
                                abort
                            not repeat:
                                record

        Returns:
            [edit df in place]
            idx_change_mapper: the mapper recording index change when repeating chain id
            is updated
        """
        idx_change_mapper = {}
        # get legal chain id
        chain_ids = set(list(map(lambda c_id: c_id.strip(), df["chain_id"])))
        if "" in chain_ids:
            chain_ids.remove("")
        legal_ids = cls._get_legal_pdb_chain_ids(chain_ids)
        ## split the df into chain dfs
        ter_line_ids = list(ter_df["line_idx"])
        chain_dfs = split_df_base_on_column_value(df, "line_idx", ter_line_ids)
        ## get a list of (loc, value) for editing in the original df
        result_loc_map = []
        recorded_chain_ids = []
        chain: pd.DataFrame
        for chain in chain_dfs:
            # check for redundant chains
            if len(chain) == 0:
                _LOGGER.debug("Found empty chain. ")
                continue
            # check existing chain ids in this chain
            existing_chain_ids = set(list(map(lambda c_id: c_id.strip(), chain["chain_id"])))
            existing_chain_ids.discard("")
            atom_missing_c_id = chain[chain.chain_id.str.strip() == ""]
            if len(atom_missing_c_id) == 0:
                # case: no missing chain id (but potentially repeating chain id)
                # divid into ATOM and HET chain
                current_chain_records = set(list(map(lambda c_id: c_id.strip(), chain["record_name"])))

                if len(existing_chain_ids) == 1:
                    # case: only 1 chain id
                    current_chain_id = list(existing_chain_ids)[0]
                    if "ATOM" in current_chain_records:
                        # case: no missing chain id; single chain id; ATOM chain
                        if current_chain_id in recorded_chain_ids:
                            _LOGGER.error("Found the same chain id in 2 different ATOM chains. Check your PDB.")
                            sys.exit(1)
                        recorded_chain_ids.append(current_chain_id)
                    else:
                        # case: no missing chain id; single chain id; HET chain
                        if current_chain_id in recorded_chain_ids:
                            new_chain_id = legal_ids.pop()
                            _LOGGER.debug(f"Found repeating chain id in a HETATM chain: assigning a new chain: {new_chain_id}")
                            result_loc_map.extend(list(zip(chain.index, [new_chain_id] * len(chain))))
                            cls._write_idx_change_when_resolve_chain_id(idx_change_mapper, chain, new_chain_id)
                        else:
                            recorded_chain_ids.append(current_chain_id)
                    continue
                else:
                    # case: more than 1 chain id (since it cannot be 0)
                    if "ATOM" in current_chain_records:
                        if not allow_multichain_in_atom:
                            # ATOM chain: multiple chain id not allowed
                            _LOGGER.error("Found multiple chain id in ATOM chain. Not allowed.")
                            sys.exit(1)
                        else:
                            # ATOM chain: (user option) multiple chain id allowed
                            _LOGGER.warning(
                                "Found multiple chain id in 1 ATOM chain. Allowed by the user option. Treating them as individual chains")
                            # case: no missing chain ids; multiple chain id; ATOM chain (user option)
                            chains_in_chain = chain.groupby("chain_id", sort=False)
                            for chain_id_in_chain, chain_in_chain in chains_in_chain:
                                if chain_id_in_chain in recorded_chain_ids:
                                    _LOGGER.error("Found the repeating chain id in a ATOM chain with multiple chain id. Check your PDB.")
                                    sys.exit(1)
                                else:
                                    recorded_chain_ids.append(chain_id_in_chain)
                            continue
                    else:
                        # HETATM chain: multiple chain id allowed and should be treated as individual chains
                        _LOGGER.warning("Found multiple chain id in 1 HETATM chain. Be careful. Treating them as individual chains")
                        # case: no missing chain ids; multiple chain id; HET chain
                        chains_in_chain = chain.groupby("chain_id", sort=False)
                        for chain_id_in_chain, chain_in_chain in chains_in_chain:
                            if chain_id_in_chain in recorded_chain_ids:
                                new_chain_id = legal_ids.pop()
                                _LOGGER.debug(f"Found repeating chain id in a HETATM chain: assigning a new chain: {new_chain_id}")
                                result_loc_map.extend(list(zip(chain_in_chain.index, [new_chain_id] * len(chain_in_chain))))
                                cls._write_idx_change_when_resolve_chain_id(idx_change_mapper, chain_in_chain, new_chain_id)
                                # no need to record since it cant be the same as any existing or new ones
                            else:
                                recorded_chain_ids.append(chain_id_in_chain)
                        continue
            else:
                if len(existing_chain_ids) == 0:
                    # case: all atoms are missing chain id
                    result_loc_map.extend(list(zip(atom_missing_c_id.index, [legal_ids.pop()] * len(atom_missing_c_id))))
                elif len(existing_chain_ids) == 1:
                    # case: only some atoms are missing chain id
                    current_chain_id = list(existing_chain_ids)[0]
                    if current_chain_id in recorded_chain_ids:
                        # divid into ATOM and HET chain
                        current_chain_records = set(list(map(lambda c_id: c_id.strip(), chain["record_name"])))
                        if "ATOM" in current_chain_records:
                            _LOGGER.error("Found the same chain id in 2 different ATOM chains. Check your PDB.")
                            sys.exit(1)
                        else:
                            new_chain_id = legal_ids.pop()
                            _LOGGER.debug(f"Found repeating chain id in a HETATM chain: assigning a new chain: {new_chain_id}")
                            result_loc_map.extend(list(zip(chain.index, [new_chain_id] * len(chain))))
                            cls._write_idx_change_when_resolve_chain_id(idx_change_mapper, chain, new_chain_id)
                    else:
                        result_loc_map.extend(list(zip(atom_missing_c_id.index, [current_chain_id] * len(atom_missing_c_id))))
                else:
                    # case: more than 1 chain id in chain
                    _LOGGER.error("Found multiple chain id together with missing chain id in 1 chain. Impossible to solve.")
                    sys.exit(1)
        # add missing chain id
        batch_edit_df_loc_value(df, result_loc_map, "chain_id")
        return idx_change_mapper

    @staticmethod
    def _get_legal_pdb_chain_ids(taken_ids) -> List[str]:
        """
        Small helper method that determines the legal PDB chain names given a list of taken ids.
        Uses all of the available 26 capitalized letters and returns in reverse order. Returns an 
        empty list when all are occupied.
        TODO in the case parsing MD results there will be 10000+ water as individual chains but still not enough for water case
             use numbers as place holders, but its making things slower, need a better solution
        Arg:
            taken_ids: iterable that contain all taken chain ids as str
        Return:
            the legal chain name list in a reversed order (for doing .pop())
        """
        result = list(string.ascii_uppercase) + list(map(lambda x: str(x), range(50000)))
        result = list(filter(lambda s: s not in taken_ids, result))
        return list(reversed(result))

    @staticmethod
    def _write_idx_change_when_resolve_chain_id(idx_change_mapper: dict, df: pd.DataFrame, new_chain_id: str):
        """specialize function for _resolve_missing_chain_id. write idx changes
        in the dataframe to a mapper"""
        for c_r_id, res in df.groupby(["chain_id", "residue_number"]):
            idx_change_mapper[c_r_id] = (new_chain_id, c_r_id[1])

    @staticmethod
    def _resolve_alt_loc(df: pd.DataFrame, keep: str = "first") -> None:
        """
        Resolves atoms with the alt_loc records
        Optional argument "keep" specifies the resolution method. Default
        is "first" which forces keeping atoms in the frist residue (earlist 
        in the file) out of multiple ones. Otherwise, keeps the specific 
        alt_loc identifier (e.g.: B), assuming it exists.

        Only one record out of multiple alt_loc is allowed. Delete rest 
        df lines in place. 
        TODO support residue specific keep strageties
        """
        alt_loc_atoms_df = df[df["alt_loc"].str.strip() != ""]
        # san check
        if len(alt_loc_atoms_df) == 0:
            _LOGGER.debug("No alt_loc to resolve.")
            return
        # solve
        # get a list of "loc" for deleting in the original df
        delete_loc_list = []
        # treat in residues
        alt_loc_residues = alt_loc_atoms_df.groupby(["residue_number", "chain_id"], sort=False)
        for r_id_c_id, res_df in alt_loc_residues:
            alt_res_dfs = res_df.groupby("alt_loc", sort=False)
            if len(alt_res_dfs) == 1:
                _LOGGER.debug(f"Only 1 alt_loc id found in residue {r_id_c_id[1], r_id_c_id[0]}. No need to resolve")
                continue
            _LOGGER.debug(f"Dealing with alt_loc residue: {r_id_c_id[1], r_id_c_id[0]}")
            if keep == "first":
                delele_res_df_locs = list(alt_res_dfs.groups.values())
                del delele_res_df_locs[0]
                for delete_lines in delele_res_df_locs:
                    delete_loc_list.extend(list(delete_lines))
            else:
                delele_res_dfs_mapper = alt_res_dfs.groups
                assert keep in delele_res_dfs_mapper
                del delele_res_dfs_mapper[keep]
                for delete_lines in list(delele_res_dfs_mapper.values()):
                    delete_loc_list.extend(list(delete_lines))

        # delete in original df
        _LOGGER.debug(f"deleting df row: {delete_loc_list}")
        df.drop(index=delete_loc_list, inplace=True)

    @staticmethod
    def _build_atoms(df: pd.DataFrame) -> Dict[str, Atom]:
        """
        build atom objects from a "standard" dataframe
        * HETATM should be in seperate chains! HETATM was meant for atoms in the small molecule
          so they cannot be in the same chain as the ATOM atom https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        Return:
            {("chain_id", "residue_idx", "residue_name", "record_name") : [Atom_obj, ...], ...}
        """
        atom_mapper = defaultdict(list)
        df.sort_values("line_idx", inplace=True)
        for i, row in df.iterrows():
            atom_obj = Atom(row)
            residue_key = (row["chain_id"].strip(), row["residue_number"], row["residue_name"].strip())
            atom_mapper[residue_key].append(atom_obj)
        return atom_mapper

    @classmethod
    def _build_residues(cls,
                        atom_mapper: Dict[str, Atom],
                        add_solvent_list: List = None,
                        add_ligand_list: List = None) -> Dict[str, Residue]:  #TODO add test for this
        """
        build Residue() objects from atom_mapper
        Arg:
            atom_mapper
            add_solvent_list: (used for categorize residues) additional names for solvent
            add_ligand_list: (used for categorize residues) additional names for ligands
            * solvent list have higher pirority
        Return:
            {"chain_id" : [Residue, ...]}
        """
        result_mapper = defaultdict(list)
        for res_key, atoms in atom_mapper.items():
            res_obj = Residue(int(res_key[1]), res_key[2], atoms)
            result_mapper[res_key[0]].append(res_obj)  # here it is {"chain_id": Residue()}
        # categorize_residue
        cls._categorize_pdb_residue(result_mapper, add_solvent_list, add_ligand_list)

        return result_mapper

    @staticmethod
    def _categorize_pdb_residue(residue_mapper: Dict[str, List[Residue]],
                                add_solvent_list: list = None,
                                add_ligand_list: list = None) -> Union[Residue, Ligand, Solvent, MetalUnit]:  #TODO add test for this
        """
        Categorize Residue base on it 3-letter name and chain info in PDB format
        Takes a mapper of {chain_id, [Residue, ...]} and converts
        them into its specialized Residue() inherited class.
        """
        if add_solvent_list is None:
            add_solvent_list = []
        if add_ligand_list is None:
            add_ligand_list = []
        for chain_id, residues in residue_mapper.items():
            peptide_chain = 0
            for residue in residues:
                # start from unknown type: canonical, non-canonical, metal, ligand, solvent
                # if canonical
                if residue.name in chem.THREE_LETTER_AA_MAPPER:
                    residue.rtype = chem.ResidueType.CANONICAL
                    peptide_chain = 1
            # non-canonical and ligand: they are only differet by if they appears
            # in a chain with canonical aa or a individual chain
            if peptide_chain:
                # only non-canonical aa can be in a peptide chain
                for residue in residues:
                    if residue.rtype == chem.ResidueType.UNKNOWN:
                        if residue.name in list(chem.METAL_MAPPER.keys()) + chem.RD_SOLVENT_LIST + add_solvent_list:
                            _LOGGER.error(
                                f"a metal or solvent residue name is found in an peptide chain {chain_id}: {residue.idx} {residue.name}")
                            sys.exit(1)
                        # TODO maybe a class for non-canonical aa in the future
                        _LOGGER.debug(f"found noncanonical {chain_id} {residue.idx}")
                        residue.rtype = chem.ResidueType.NONCANONICAL
                continue
            # non-peptide chain
            for i, residue in enumerate(residues):
                # if metal
                if residue.name in chem.METAL_MAPPER:
                    _LOGGER.debug(f"found metal {chain_id} {residue.idx}")
                    residue_mapper[chain_id][i] = residue_to_metal(residue)
                    continue
                # if solvent
                if residue.name in chem.RD_SOLVENT_LIST + add_solvent_list:
                    _LOGGER.debug(f"found solvent {chain_id} {residue.idx}")
                    residue_mapper[chain_id][i] = residue_to_solvent(residue)
                    continue
                if residue.name in (set(chem.RD_NON_LIGAND_LIST) - set(add_ligand_list)):
                    _LOGGER.debug(f"found trash {chain_id} {residue.idx}")
                    residue.rtype = chem.ResidueType.TRASH
                    continue
                # the left cases are ligands
                _LOGGER.debug(f"found ligand {chain_id} {residue.idx}")
                residue_mapper[chain_id][i] = residue_to_ligand(residue)

    @staticmethod
    def _build_chains(res_mapper: Dict[str, Residue], remove_trash: bool) -> Dict[str, Chain]:
        """
        build Chain() objects from residue mapper
        """
        chain_list = []
        for chain_id, residues in res_mapper.items():
            ch_obj = Chain(chain_id, residues)  # TODO decouple chain name with pdb chain id? How to solve too many water chain
            if remove_trash:
                ch_obj.remove_trash()
            chain_list.append(ch_obj)
        return chain_list

    #endregion

    @classmethod
    def _write_pdb_chain(cls, chain: Chain) -> str:
        """
        make the file string for a pdb chain record
        """
        result = ""
        chain.sort_residues()
        for res in chain:
            res: Residue
            result += cls._write_pdb_residue(res)
        result += f"TER{os.linesep}"
        return result

    @classmethod
    def _write_pdb_residue(cls, res: Residue) -> str:
        """
        make the file string for a pdb residue record
        """
        result = ""
        res.sort_atoms()
        for atom in res:
            atom: Atom
            result += cls._write_pdb_atom(atom)
        return result

    @staticmethod
    def _write_pdb_atom(atom: Atom) -> str:
        """
        make the file string for a pdb ATOM record
        the function intepret information from
        https://www.wwpdb.org/documentation/file-format
        """
        #build an amber style line here
        l_type = f"{'ATOM':<6}"
        a_index = f"{atom.idx:>5d}"

        if len(atom.name) > 3:
            a_name = f"{atom.name:<4}"
        else:
            a_name = f"{atom.name:<3}"
            a_name = " " + a_name
        r_name = f"{atom.residue.name:>3}"

        c_index = atom.residue.chain.name
        r_index = f"{atom.residue.idx:>4d}"
        x = f"{atom.coord[0]:>8.3f}"
        y = f"{atom.coord[1]:>8.3f}"
        z = f"{atom.coord[2]:>8.3f}"

        alt_loc_id = f"{'':1}"
        insert_code = f"{'':1}"
        occupancy = f"{1.00:>6.2f}"

        if atom.b_factor is None:
            temp_factor = f"{0.0:>6.2f}"
        else:
            temp_factor = f"{atom.b_factor:>6.2f}"

        seg_id = f"{'':<4}"

        if atom.element is None:
            element = f"{'':>2}"
        else:
            element = f"{atom.element:>2}"

        if atom.charge is None:
            charge = f"{'':2}"
        else:
            charge = f"{atom.charge:2}"

        #example: ATOM   5350  HB2 PRO   347      32.611  15.301  24.034  1.00  0.00
        line = f"{l_type}{a_index} {a_name}{alt_loc_id}{r_name} {c_index}{r_index}{insert_code}   {x}{y}{z}{occupancy}{temp_factor}      {seg_id}{element}{charge}{os.linesep}"
        return line

    @dispatch
    def _(self):
        """
        dummy method for dispatch
        """
        pass


#region == util of PDB file handling ==
# these tools handles those tasks that dont worth converting to Structure


# index related
def restore_pdb_index(index_mapper: Dict[tuple, tuple], target_pdb: str, out_path: str) -> None:
    """Helper method that renumbers the residue chain id's and residue numbers in the
    new structure found in target_pdb according to the index mapper (go back to old index)
    Args:
        index_mapper:
            The mapper for replacing the index in target pdb
            {(new_chain_id, new_residue_idx):(old_chain_id, old_residue_idx),
            ...}
        target_pdb: The path of the target .pdb file.
        out_path: the path for the output pdb
    """
    pdb_panda_obj = PandasPdb().read_pdb(target_pdb)
    pdb_df = get_pdb_df(target_pdb)  # also fix for chain id
    for index, row in pdb_df.iterrows():
        old_chain_id, old_resi_idx = index_mapper[(row["chain_id"], row["residue_number"])]
        # edit index in place
        pdb_df.loc[index, "chain_id"] = old_chain_id
        pdb_df.loc[index, "residue_number"] = old_resi_idx
    pdb_panda_obj.df["ATOM"] = pdb_df
    pdb_panda_obj.to_pdb(out_path, records=["ATOM", "OTHERS"])


def get_index_mapper_from_pdb(npdb: str, opdb: str, method: Union[None, str] = None, **kwargs) -> Dict[tuple, tuple]:
    """get the index mapper from matching old pdb and new pdb
    (new pdb -> old pdb)"""
    npdb_df = get_pdb_df(npdb)
    opdb_df = get_pdb_df(opdb)
    if method is None:
        method = "default"
    result = INDEX_MAPPING_METHODS[method](npdb_df, opdb_df, **kwargs)
    return result


def pdb_index_mapping_by_order(npdb_df: pd.DataFrame, opdb_df: pd.DataFrame) -> Dict[tuple, tuple]:
    """default method of mapping indexes (new pdb -> new pdb)
    Just align the indexes by orders in the PDB file:
        The residue of the same sequence appearing in the PDB file should be mapped
        as the same residue and thus using the correponding opdb index.

    WARNING: those do not consider the case that orders of residue sequence in the file
    can be inverted after some process. (like pymol or enzyhtp they will invert the
    order base on residue indexes)"""
    nkeys, okeys = get_pdb_index_key(npdb_df), get_pdb_index_key(opdb_df)
    if len(okeys) != len(nkeys):
        _LOGGER.error(f"Residue key number does not match between old and new pdb. Index mapping is impossible")
        raise IndexMappingError
    if if_list_contain_repeating_element(okeys):
        _LOGGER.error(f"Found repeating residue key in old pdb. Index mapping is impossible")
        raise IndexMappingError
    if if_list_contain_repeating_element(nkeys):
        _LOGGER.error(f"Found repeating residue key in new pdb. Index mapping is impossible")
        raise IndexMappingError

    idx_mapper = dict(zip(nkeys, okeys))
    return idx_mapper


@dispatch
def get_pdb_index_key(pdb: str) -> List[Tuple[str, str]]:
    """get a index key of every residue in the pdb in the form of
    (chain_id, residue_number) with the same order in the PDB file"""
    pdb_df = get_pdb_df(pdb)
    result = list_remove_adjacent_duplicates(list(zip(pdb_df.residue_name, pdb_df.chain_id, pdb_df.residue_number)))
    result = [(i, j) for _, i, j in result]
    return result


@dispatch
def get_pdb_index_key(pdb_df: pd.DataFrame) -> List[Tuple[str, str]]:
    """get a index key of every residue in the pdb in the form of
    (chain_id, residue_number) with the same order in the PDB file"""
    result = list_remove_adjacent_duplicates(list(zip(pdb_df.residue_name, pdb_df.chain_id, pdb_df.residue_number)))
    result = [(i, j) for _, i, j in result]
    return result


INDEX_MAPPING_METHODS = {
    "default": pdb_index_mapping_by_order,
    "by_order": pdb_index_mapping_by_order,
}


# general
def get_pdb_df(pdb: str, fix_chain_id: bool = True) -> pd.DataFrame:
    """get the df of pdb. Fix chain id by default.
    The idea is chain id information is already implicitly defined by TER lines
    Fixing it is just make it explicit"""
    df: pd.DataFrame = PandasPdb().read_pdb(pdb).df
    result = pd.concat((df["ATOM"], df["HETATM"]), ignore_index=True)
    result.sort_values("line_idx", inplace=True)
    if fix_chain_id:
        pdb_ter_df = df["OTHERS"][df["OTHERS"].record_name == "TER"]
        PDBParser._resolve_missing_chain_id(result, pdb_ter_df)
    return result


#endregion
