"""Generation/construction of Structure objects from mmCIF files and exporting it vice versa
Definition of the mmCIF file format (https://www.wwpdb.org/documentation/file-format-content/cif33/cc.html)
Also contains util function that requires knoweledge of the mmCIF format

`Structure` enrichment wishlist (to keep the mmCIF parsing advantages while staying
PDB compatible. the parser consumes the richer mmCIF information at the parse time
and the `Structure` data model is the bottleneck of how much of it survives):
    1. track label_entity_id/label_seq_id of each residue (e.g.: as optional
       `Residue` fields). Currently the polymer-membership signal is consumed
       once to categorize residues and then dropped, so the exported file
       hardcodes entity '1' and label_seq '.' (see _write_cif_atom_site_row()).
    2. track the original ATOM/HETATM record of each residue. Currently the
       writer re-derives the record from the residue category (canonical or
       modified -> ATOM), which normalizes unusual files instead of faithfully
       round-tripping them.
    3. member views for the mixed-content chains (mmCIF keeps waters/metals/
       ligands on the auth chain id of the polymer): rtype-filtered residue
       views (e.g.: polymer-only) and a polymer-only Chain.sequence, so the
       downstream modules do not see 'HOH' in a sequence string.
    4. resolve the (chain_id, resi) key collision in
       Structure.find_residue_with_key() explicitly (e.g.: warn) instead of
       silently taking the first match. Waters renumbered from 1 on the same
       auth chain would shadow the polymer residue with the same number.
    5. multi-model support: Structure` only holds the selected model.

Author: Zhong, Yinjie <yinjie.zhong@outlook.com>
Date: 2026-10-07
"""
import copy
import os
from collections import defaultdict
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from plum import dispatch
from biopandas.mmcif import PandasMmcif

from enzy_htp.core import _LOGGER
import enzy_htp.core.file_system as fs
import enzy_htp.chemical as chem
from ._interface import StructureParserInterface
from .pdb_io import PDBParser
from ..atom import Atom
from ..metal_atom import residue_to_metal
from ..residue import Residue
from ..solvent import residue_to_solvent
from ..ligand import residue_to_ligand
from ..modified_residue import residue_to_modified_residue
from ..chain import Chain
from ..structure import Structure


class CIFParser(StructureParserInterface):
    """
    Parser for covert mmCIF into Structure and vice versa
    * Different from PDB, mmCIF allows non-polymer (HETATM) components to share
      the same chain id (auth_asym_id) with the polymer (ATOM) ones (e.g.: water
      with the auth chain id of the protein). They are merged into the same `Chain` instance.
    * Consistency with PDBParser is verified on the test data pairs (see
      `test/structure/structure_io/test_cif_io.py`):
        - 3FCR: all 3829 atoms are recovered from the mmCIF with identical
          coordinates and the same categorization (0 mismatch over 777 residues).
          The modified residues (MSE x6, LLP x1) written as HETATM records
          sharing the polymer entity id are recognized via the
          label_entity_id/label_seq_id signals.
        - 1NVG: all 2616 atoms of the legacy de-watered PDB file are recovered
          from the mmCIF with identical coordinates, b-factors, and elements
          and the same categorization (0 mismatch over 349 residues, incl. the
          2 ZN ions -> MetalUnit). The mmCIF side keeps strictly more
          information than its PDB counterpart: the 6 alt loc atoms of HIS 39
          keep the kept alternates' alt loc 'A' designation and the real
          occupancy (0.6) where the PDB file has both stripped (alternates
          removed, occupancy flattened to 1.0); the 120 crystallographic waters
          are kept as `Solvent` instance on the auth chain A of the protein where the
          PDB file has none; and nothing is exiled to a separate chain (the
          PDB file puts the ZN ions on their own chain B) since mmCIF carries
          the auth chain id on every atom.
    """

    def __init__(self) -> None:  # pylint: disable=super-init-not-called
        """
        not sure what can be init here. Maybe some configuration?
        """
        pass

    # knowledge
    ## the tags of the _atom_site loop in the order used by the RCSB.
    ## This is also the read/write interface of this parser.
    ATOM_SITE_TAGS = [
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
        "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
        "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
        "B_iso_or_equiv", "pdbx_formal_charge", "auth_seq_id", "auth_comp_id",
        "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"
    ]
    ## the data block name of the exported file since `Structure` does not track an entry id
    DATA_BLOCK_NAME = "enzy_htp"

    # interface
    @classmethod
    def get_structure(cls,
                      path: str,
                      model: int = 0,
                      add_solvent_list: List = None,
                      add_ligand_list: List = None,
                      remove_trash: bool = True,
                      give_idx_map: bool = False,
                      alt_loc_keep: str = "first",
                      single_alt_loc: str = "clean",) -> Union[Structure, tuple]:
        """
        Converting a mmCIF file (as its path) into the a `Structure` instance.

        Arg:
            path:
                the file path of the mmCIF file
            model:
                The selected model index if there are multiple models in the file
                (start from 0 default: 0)
                (models are resolved by the "pdbx_PDB_model_num" tag)
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
                if return a tuple of (Structure, idx_change_mapper). For mmCIF the
                chain id (auth_asym_id) is carried by every atom natively, no chain
                id resolution ever happens and thus the mapper is always an empty
                dict. It is kept here for interface consistency with PDBParser.
            alt_loc_keep:
                specifies the resolution method for alt_loc records.
                Default is "first" which keeps atoms with the lexicographically smallest
                alt_loc identifier (e.g., A before B) for each residue, independent of
                DataFrame row order. Otherwise, if a specific alt_loc identifier is
                provided (e.g., "B"), it keeps that identifier, assuming it exists. The "else"
                operator can be used to specify a fallback alt_loc identifier in cases that
                the preferred one does not exist in all residues that contains alt_loc
                (e.g., "B else A" keeps "B" when it exists in a residue, otherwise keeps "A").
            single_alt_loc:
                specifies how to treat the case when there is only one alt_loc in a residue.
                Default is "clean", which deletes the alt_loc identifier (sets it to blank)
                and keeps the atoms. Other options are "keep" (keeps the alt_loc identifier
                as is) and "delete" (deletes the atoms with alt_loc identifier).

        Return:
            Structure
        """
        _LOGGER.debug(f"working on {path}")
        cls._check_valid_cif(path)
        # covert to dataframe with biopanda (maybe write our own in the future)
        input_cif = PandasMmcif()
        input_cif.read_mmcif(path)

        #region (mmCIF conundrums)
        # deal with multiple model
        df_parts = [input_cif.df[key] for key in ("ATOM", "HETATM") if key in input_cif.df and len(input_cif.df[key]) > 0]
        target_model_df = cls._get_target_model(pd.concat(df_parts, ignore_index=True), model=model)
        # normalize the df into the biopandas pdb naming so the downstream
        # building steps can be shared with PDBParser
        target_model_df = cls._normalize_atom_site_df(target_model_df)
        cls._resolve_alt_loc(target_model_df, keep=alt_loc_keep, single_alt_loc=single_alt_loc)  # resolve alt loc record in place by delele redundant df rows
        # Unlike PDB, mmCIF carries the chain id (auth_asym_id) of every atom
        # natively. No missing/repeating chain id resolution is needed.
        #endregion (Add here to address more problem related to the mmCIF file format here)

        # target_model_df below should be "standard" --> start building
        # reuse the format agnostic parts of the PDB building steps.
        # see CIFParser._normalize_atom_site_df() for the normalization making this possible
        atom_mapper: Dict[tuple, Atom]
        res_info_mapper: Dict[tuple, tuple]
        atom_mapper, res_info_mapper = cls._build_atoms(target_model_df)
        res_mapper = cls._build_residues(atom_mapper, res_info_mapper, add_solvent_list, add_ligand_list)
        chain_list: Dict[str, Chain] = PDBParser._build_chains(res_mapper, remove_trash)

        if give_idx_map:
            return Structure(chain_list), {}
        return Structure(chain_list)

    @classmethod
    def save_structure(cls,
                       outfile: str,
                       stru: Structure,
                       if_renumber: bool = True,
                       if_fix_atomname: bool = True,) -> str:
        """Inverse of CIFParser.get_structure(). Given a `Structure` instance, save it to the given .cif path.

        Args:
            outfile: Path to save the `Structure` instance as a string.
            stru: The `Structure` instance to save.
            if_renumber: Should atoms be renumbered from 1?
            if_fix_atomname: Should atoms be renamed to PDB convention?

        Returns:
            Path to the saved `Structure` instance as a string.
        """
        if not isinstance(outfile, str):
            raise TypeError("outfile must be a string representing the file path.")

        content: str = cls.get_file_str(stru, if_renumber, if_fix_atomname)

        fs.write_lines(outfile, content.splitlines())

        return outfile

    @classmethod
    @dispatch
    def get_file_str(
            cls,
            stru: Chain,
            if_renumber: bool = True,
            if_fix_atomname: bool = True,
        ) -> str:  # pylint: disable=function-redefined
        """
        dispatch for supporting get cif file str with Chain only
        """
        if if_renumber:
            stru.parent.renumber_atoms()
        result_str = cls._write_data_block_header() + cls._write_atom_site_loop_header()
        result_str += cls._write_cif_chain(stru)
        result_str += f"#{os.linesep}"
        return result_str

    @classmethod
    @dispatch
    def get_file_str(
            cls,
            stru: Residue,
            if_renumber: bool = True,
            if_fix_atomname: bool = True,
            ) -> str:  # pylint: disable=function-redefined
        """
        dispatch for supporting get cif file str with Residue only
        """
        if if_renumber:
            stru.parent.parent.renumber_atoms()
        result_str = cls._write_data_block_header() + cls._write_atom_site_loop_header()
        result_str += cls._write_cif_residue(stru)
        result_str += f"#{os.linesep}"
        return result_str

    @classmethod
    @dispatch
    def get_file_str(
            cls,
            stru: Atom,
            if_renumber: bool = True,
            if_fix_atomname: bool = True,
        ) -> str:  # pylint: disable=function-redefined
        """
        dispatch for supporting get cif file str with Atom only
        """
        result_str = cls._write_data_block_header() + cls._write_atom_site_loop_header()
        result_str += cls._write_cif_atom_site_row(stru)
        result_str += f"#{os.linesep}"
        return result_str

    @classmethod
    @dispatch
    def get_file_str(   # pylint: disable=function-redefined
            cls,
            stru: Structure,
            if_renumber: bool = True,
            if_fix_atomname: bool = True,
            ) -> str:
        """
        Convert `Structure` instance into mmCIF file string. Only the simplest function is need for
        enzyme modeling.

        Args:
            stru: The `Structure` instance to convert.
            if_renumber: Whether to renumber atoms starting from 1.
            if_fix_atomname: Whether to fix atom names to PDB conventions.

        TODO support fixing all atom names before writing
        TODO the file only contains the _atom_site loop. Add more categories
             (entity, struct_asym ...) when needed.
             A better way to do this is to make derivate Parser classes and overwrite this as needed.
        Return:
            a string of the mmCIF file
        """
        stru = copy.deepcopy(stru)
        stru.sort_chains()
        if if_renumber:
            stru.renumber_atoms()
        result_str = cls._write_data_block_header() + cls._write_atom_site_loop_header()
        for chain in stru:
            chain: Chain
            result_str += cls._write_cif_chain(chain)
        result_str += f"#{os.linesep}"
        return result_str

    #region == cif -> Stru ==
    @staticmethod
    def _check_valid_cif(cif_path: str) -> None:
        """Helper function that ensures the supplied cif_path contain a valid mmCIF file.
        Private to structure.structure_io.cif_io.py. Should NOT be called externally."""
        # check for right extension
        ext: str = fs.get_file_ext(cif_path)
        if ext.lower() not in (".cif", ".mmcif"):
            message = f"Supplied file '{cif_path}' is NOT a mmCIF file."
            _LOGGER.error(message)
            raise ValueError(message)
        #check for file existance
        if not os.path.exists(cif_path):
            message = f"Supplied file '{cif_path}' does NOT exist."
            _LOGGER.error(message)
            raise ValueError(message)
        # check for character encoding
        for idx, ll in enumerate(fs.lines_from_file(cif_path)):
            if not ll.isascii():
                message = f"The mmCIF '{cif_path}' contains non-ASCII text and is invalid in line {idx}: '{ll}'."
                _LOGGER.error(message)
                raise ValueError(message)

    @staticmethod
    def _get_target_model(df: pd.DataFrame, model: int) -> pd.DataFrame:
        """
        get the #model model part of the dataframe when there are multiple models
        in the mmCIF file.
        Arg:
            df: the combined ATOM/HETATM dataframe of a mmCIF file
            model: the index of target model
                (start from 0)
                (models are resolved by the "pdbx_PDB_model_num" tag)
        Return: a copy of Atom and HETATM records from the input
                dataframe correponding to the target model
        """
        model_ids = sorted(df["pdbx_PDB_model_num"].unique())
        if model < 0 or model >= len(model_ids):
            message = f"Model index {model} is out of range. The file only has {len(model_ids)} model(s)."
            _LOGGER.error(message)
            raise ValueError(message)
        return df[df["pdbx_PDB_model_num"] == model_ids[model]].copy()

    @staticmethod
    def _normalize_atom_site_df(df: pd.DataFrame) -> pd.DataFrame:
        """
        normalize the biopandas mmcif atom_site dataframe into the biopandas pdb naming
        so that the downstream building steps can be shared with PDBParser
        (Atom.from_biopandas(), PDBParser._build_atoms(), PDBParser._resolve_alt_loc() ...)
        * author provided ids (auth_*) are preferred over the label_* ones since
          they align with the PDB format convention that EnzyHTP uses.
        * the placeholder value of the CIF format ('.', '?') and the 'None' string
          converted by biopandas are all mapped to the blank convention of the
          biopandas pdb dataframe (empty string or NaN)
        """
        def pick_column(preferred: str, fallback: str = None) -> pd.Series:
            """get the column of the preferred tag or use the fallback one"""
            if preferred in df.columns:
                return df[preferred]
            if fallback is not None and fallback in df.columns:
                _LOGGER.debug(f"mmCIF atom_site is missing the {preferred} tag. Using {fallback} instead.")
                return df[fallback]
            message = f"mmCIF atom_site is missing both the {preferred} and the {fallback} tag."
            _LOGGER.error(message)
            raise ValueError(message)

        def clean_str(series: pd.Series) -> pd.Series:
            """map the cif placeholder value ('.', '?') and the 'None' string converted
            by biopandas into the empty string"""
            return series.fillna("").astype(str).replace({".": "", "?": "", "None": ""})

        def clean_number(series: pd.Series) -> pd.Series:
            """map the cif placeholder value ('.', '?') and the 'None' string converted
            by biopandas into NaN"""
            return pd.to_numeric(clean_str(series).replace({"": None}), errors="coerce")

        def parse_charge(series: pd.Series) -> pd.Series:
            """convert the formal charge notation of CIF (e.g.: '1-', '-2') into float"""
            def convert(value: str) -> float:
                if value in ('.', '?', '', 'None', None):
                    return np.nan
                try:
                    if value.endswith("-"):
                        return -float(value[:-1])
                    if value.endswith("+"):
                        return float(value[:-1])
                    return float(value)
                except ValueError:
                    _LOGGER.debug(f"Unparsable formal charge '{value}'. Setting it to NaN.")
                    return np.nan
            return series.map(convert)

        result = pd.DataFrame()
        result["record_name"] = df["group_PDB"]
        result["atom_number"] = clean_number(df["id"])
        result["atom_name"] = clean_str(pick_column("auth_atom_id", "label_atom_id"))
        result["residue_name"] = clean_str(pick_column("auth_comp_id", "label_comp_id")).str.strip()
        result["chain_id"] = clean_str(pick_column("auth_asym_id", "label_asym_id")).str.strip()
        result["residue_number"] = clean_number(pick_column("auth_seq_id", "label_seq_id"))
        result["x_coord"] = clean_number(df["Cartn_x"])
        result["y_coord"] = clean_number(df["Cartn_y"])
        result["z_coord"] = clean_number(df["Cartn_z"])
        result["occupancy"] = clean_number(df["occupancy"])
        result["b_factor"] = clean_number(df["B_iso_or_equiv"])
        result["alt_loc"] = clean_str(df["label_alt_id"])
        result["element_symbol"] = clean_str(df["type_symbol"])
        result["charge"] = parse_charge(clean_str(df["pdbx_formal_charge"]))
        result["line_idx"] = range(len(df))
        # carry the polymer-membership signals for the residue categorization:
        # a HETATM residue sharing the polymer entity id with a numeric
        # label_seq_id is a modified monomer in the polymer sequence.
        if "label_entity_id" in df.columns:
            result["entity_id"] = clean_str(df["label_entity_id"])
        if "label_seq_id" in df.columns:
            result["label_seq_id"] = clean_str(df["label_seq_id"])
        return result

    @classmethod
    def _resolve_alt_loc(cls, df: pd.DataFrame, keep: str = "first", single_alt_loc: str = "clean") -> None:
        """resolve alt loc records of the normalized atom_site dataframe.
        The logic is the same as the PDB format and thus reusing the PDB implementation"""
        PDBParser._resolve_alt_loc(df, keep=keep, single_alt_loc=single_alt_loc)

    @staticmethod
    def _build_atoms(df: pd.DataFrame) -> Tuple[Dict[tuple, Atom], Dict[tuple, tuple]]:
        """
        build atom objects from a "standard" dataframe. Same as the PDB one except
        the group_PDB record is kept in the residue key for the categorization step.
        (in mmCIF, polymer ATOM and non-polymer HETATM components can share the same
        chain id, which is the main difference from the PDB format)
        Return:
            atom_mapper:
                {("chain_id", "residue_idx", "residue_name", "record_name") : [Atom_obj, ...], ...}
            res_info_mapper:
                {("chain_id", "residue_idx", "residue_name") : (entity_id, label_seq_id), ...}
                carrying the polymer-membership signals of each residue. ("", "") when
                the file does not provide the label_entity_id/label_seq_id tags.
        """
        atom_mapper = defaultdict(list)
        res_info_mapper = {}
        df.sort_values("line_idx", inplace=True)
        for i, row in df.iterrows():
            atom_obj = Atom.from_biopandas(row)
            residue_key = (row["chain_id"].strip(), row["residue_number"], row["residue_name"].strip(), row["record_name"].strip())
            atom_mapper[residue_key].append(atom_obj)
            res_key = residue_key[:3]
            if res_key not in res_info_mapper:
                res_info_mapper[res_key] = (row.get("entity_id", ""), row.get("label_seq_id", ""))
        return atom_mapper, res_info_mapper

    @classmethod
    def _build_residues(cls,
                        atom_mapper: Dict[tuple, Atom],
                        res_info_mapper: Dict[tuple, tuple],
                        add_solvent_list: List = None,
                        add_ligand_list: List = None) -> Dict[str, Residue]:
        """
        build `Residue` instances from atom_mapper
        Arg:
            atom_mapper
            res_info_mapper: the polymer-membership signals of each residue
                (see CIFParser._build_atoms() for the definition)
            add_solvent_list: (used for categorize residues) additional names for solvent
            add_ligand_list: (used for categorize residues) additional names for ligands
            * solvent list have higher pirority
        Return:
            {"chain_id" : [Residue, ...]}
        """
        result_mapper = defaultdict(list)
        record_mapper = {}
        for res_key, atoms in atom_mapper.items():
            res_obj = Residue(int(res_key[1]), res_key[2], atoms)
            result_mapper[res_key[0]].append(res_obj)  # here it is {"chain_id": Residue()}
            record_mapper[res_key[:3]] = res_key[3]
        # categorize_residue
        cls._categorize_cif_residue(result_mapper, record_mapper, res_info_mapper, add_solvent_list, add_ligand_list)

        return result_mapper

    @staticmethod
    def _categorize_cif_residue(residue_mapper: Dict[str, List[Residue]],
                                record_mapper: Dict[tuple, str],
                                res_info_mapper: Dict[tuple, tuple],
                                add_solvent_list: list = None,
                                add_ligand_list: list = None) -> None:
        """
        Categorize Residue base on its name, chain info, and the ATOM/HETATM record
        in mmCIF format. Takes a mapper of {chain_id, [Residue, ...]}, a mapper of
        {(chain_id, residue_idx, residue_name): record_name}, and a mapper of the
        polymer-membership signals (see CIFParser._build_atoms()) and converts
        them into its specialized `Residue` inherited classes.
        * different from PDB, mmCIF allows non-polymer (HETATM) components to share
          the same chain id with the polymer (ATOM) ones. So the categorization is
          done per residue with its own record instead of per chain.
        * a HETATM residue can still be a polymer monomer when it shares the
          polymer entity id and has a numeric label_seq_id. This is how modified
          residues appear in mmCIF files converted from the PDB format.
        """
        if add_solvent_list is None:
            add_solvent_list = []
        if add_ligand_list is None:
            add_ligand_list = []
        # collect the entity ids of the polymer (ATOM) records
        polymer_entity_ids = set()
        for res_key, record_name in record_mapper.items():
            if record_name == "ATOM":
                entity_id = res_info_mapper.get(res_key, ("", ""))[0]
                if entity_id != "":
                    polymer_entity_ids.add(entity_id)
        for chain_id, residues in residue_mapper.items():
            peptide_chain = 0
            for residue in residues:
                # if canonical
                if residue.name in chem.THREE_LETTER_AA_MAPPER:
                    residue.rtype = chem.ResidueType.CANONICAL
                    peptide_chain = 1
            for i, residue in enumerate(residues):
                if residue.rtype == chem.ResidueType.CANONICAL:
                    continue
                res_key = (chain_id, residue.idx, residue.name)
                record_name = record_mapper[res_key]
                entity_id, label_seq_id = res_info_mapper.get(res_key, ("", ""))
                if peptide_chain and (record_name == "ATOM" or (entity_id in polymer_entity_ids and label_seq_id.isdigit())):
                    # only modified aa can be a polymer residue in a peptide chain
                    # (either written as an ATOM record or as a HETATM one sharing
                    # the polymer entity with a numeric label_seq_id)
                    _LOGGER.debug(f"found modified residue {chain_id} {residue.idx}")
                    residue_mapper[chain_id][i] = residue_to_modified_residue(residue)
                    continue
                # the left cases are all non-polymer components and are categorized by
                # their names. (an ATOM-record residue in a non-peptide chain, e.g.:
                # nucleic acid, is treated like a ligand as in the PDB flow)
                # if metal
                if residue.name in chem.METAL_MAPPER:
                    if residue.name == "U":
                        _LOGGER.warning(f"found a residue named U ({chain_id} {residue.idx})."
                                        " This could lead to a bug that it may be a metal or"
                                        " a nucleotide.")
                        if residue.num_atoms > 1:
                            # make it a ligand for now. TODO: add support for NA
                            _LOGGER.warning(f"found more than 1 atom in U ({chain_id} {residue.idx})."
                                            " treating it as nucleotide. make sure this is what you want!")
                            _LOGGER.debug(f"found ligand {chain_id} {residue.idx}")
                            residue_mapper[chain_id][i] = residue_to_ligand(residue)
                            continue
                        else:
                            _LOGGER.warning(f"found only 1 atom in U ({chain_id} {residue.idx})."
                                            " treating it as metal. make sure this is what you want!")
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

    #endregion

    #region == Stru -> cif ==
    @staticmethod
    def _write_data_block_header() -> str:
        """make the file string for the data block header.
        a CIF file must start with a data block and the _entry.id category
        is required by some parser (e.g.: biopandas)"""
        result = f"data_{CIFParser.DATA_BLOCK_NAME}{os.linesep}"
        result += f"#{os.linesep}"
        result += f"_entry.id{os.linesep}"
        result += f"{CIFParser.DATA_BLOCK_NAME}{os.linesep}"
        result += f"#{os.linesep}"
        return result

    @staticmethod
    def _write_atom_site_loop_header() -> str:
        """make the file string for the _atom_site loop header"""
        result = f"loop_{os.linesep}"
        for tag in CIFParser.ATOM_SITE_TAGS:
            result += f"_atom_site.{tag}{os.linesep}"
        return result

    @classmethod
    def _write_cif_chain(cls, chain: Chain) -> str:
        """make the atom_site rows for a chain record"""
        result = ""
        chain.sort_residues()
        for res in chain:
            res: Residue
            result += cls._write_cif_residue(res)
        return result

    @classmethod
    def _write_cif_residue(cls, res: Residue) -> str:
        """make the atom_site rows for a residue record"""
        result = ""
        res.sort_atoms()
        for atom in res:
            atom: Atom
            result += cls._write_cif_atom_site_row(atom)
        return result

    @classmethod
    def _write_cif_atom_site_row(cls, atom: Atom) -> str:
        """make the atom_site row for an atom record
        the function intepret information from
        https://www.wwpdb.org/documentation/file-format-content/cif33/cc.html
        Simplified conventions of EnzyHTP:
            1. label_* ids are the same as the auth_* ones except label_seq_id and
               label_entity_id which are not tracked by `Structure` ('1' and '.')
            2. the insertion code is not tracked by `Atom` and always written as '.'
            3. all atoms are written under model #1
            4. residues that are not canonical or modified are written as HETATM
        """
        res: Residue = atom.residue
        chain: Chain = res.chain
        group = "ATOM" if (res.is_canonical() or res.is_modified()) else "HETATM"
        occupancy = f"{atom.occupancy:.2f}" if atom.occupancy is not None else "1.00"
        b_factor = f"{atom.b_factor:.2f}" if atom.b_factor is not None else "0.00"
        if atom.charge is not None and int(atom.charge) != 0:
            charge = f"{int(atom.charge):>2}"
        else:
            charge = f"{'?':>2}"
        alt_loc = atom.alt_loc if atom.alt_loc else "."
        fields = [
            f"{group:<6}",                                  # group_PDB
            f"{atom.idx:>5d}",                              # id
            f"{atom.element:>2}",                           # type_symbol
            f"{cls._cif_quote(atom.name):<6}",              # label_atom_id
            f"{alt_loc}  ",                                 # label_alt_id
            f"{cls._cif_quote(res.name):<5}",               # label_comp_id
            f"{cls._cif_quote(chain.name):<3}",             # label_asym_id
            "1    ",                                        # label_entity_id
            ".    ",                                        # label_seq_id
            ".   ",                                         # pdbx_PDB_ins_code
            f"{atom.coord[0]:>8.3f}",                       # Cartn_x
            f"{atom.coord[1]:>8.3f}",                       # Cartn_y
            f"{atom.coord[2]:>8.3f}",                       # Cartn_z
            f"{occupancy:>6}",                              # occupancy
            f"{b_factor:>6}",                               # B_iso_or_equiv
            f"{charge}  ",                                  # pdbx_formal_charge
            f"{res.idx:>5d}",                               # auth_seq_id
            f"{cls._cif_quote(res.name):<5}",               # auth_comp_id
            f"{cls._cif_quote(chain.name):<3}",             # auth_asym_id
            f"{cls._cif_quote(atom.name):<6}",              # auth_atom_id
            "1",                                            # pdbx_PDB_model_num
        ]
        return " ".join(fields) + os.linesep

    @staticmethod
    def _cif_quote(value: str) -> str:
        """quote a value into a CIF string when it contains characters not allowed
        in a bare CIF token (whitespace, quotes, '#', '$', ';', '[', ']').
        An empty value is quoted into '' to stay a valid token."""
        value = str(value)
        if value != "" and not any(c in value for c in " \t#'\"$_;[]"):
            return value
        if '"' not in value:
            return f'"{value}"'
        if "'" not in value:
            return f"'{value}'"
        # last resort: the value contains both quote types. CIF does not support
        # escaping quotes in a quoted string. replace double quotes with single ones
        _LOGGER.warning(f"CIF value '{value}' contains both quote types. Replacing double quotes with single ones.")
        return f'"{value.replace(chr(34), chr(39))}"'

    #endregion

    @dispatch
    def _(self):
        """
        dummy method for dispatch
        """
        pass


#region == util of mmCIF file handling ==
# these tools handles those tasks that dont worth converting to Structure


def get_cif_df(cif: str, model: int = 0) -> pd.DataFrame:
    """get the normalized df of the atom_site loop of a mmCIF file.
    The df follows the biopandas pdb naming (see CIFParser._normalize_atom_site_df())
    Arg:
        cif: the file path of the mmCIF file
        model: the index of target model (start from 0)
    """
    cif_obj = PandasMmcif()
    cif_obj.read_mmcif(cif)
    df_parts = [cif_obj.df[key] for key in ("ATOM", "HETATM") if key in cif_obj.df and len(cif_obj.df[key]) > 0]
    target_model_df = CIFParser._get_target_model(pd.concat(df_parts, ignore_index=True), model=model)
    return CIFParser._normalize_atom_site_df(target_model_df)


#endregion
