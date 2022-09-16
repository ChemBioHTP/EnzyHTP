"""Generation/construction of Structure objects from PDB files and exporting it vice versa
Definition of PDB file format (http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)

Author: shaoqz, <shaoqz@icloud.com>
Date: 2022-08-01
"""
from collections import defaultdict
import os
import string
import sys
from typing import Dict, List, Tuple, Union
from biopandas.pdb import PandasPdb
import pandas as pd


from ._interface import StructureParserInterface
from ..atom import Atom
from ..metal_atom import MetalUnit, residue_to_metal
from ..residue import Residue
from ..solvent import Solvent, residue_to_solvent
from ..ligand import Ligand, residue_to_ligand
from ..chain import Chain
from ..structure import Structure
from enzy_htp.core import _LOGGER
import enzy_htp.core.file_system as fs
import enzy_htp.chemical as chem

# pylint: disable=logging-fstring-interpolation
class PDBParser(StructureParserInterface):
    '''
    Parser for covert PDB into Structure and vice versa
    '''
    def __init__(self) -> None: # pylint: disable=super-init-not-called
        '''
        not sure what can be init here. Maybe some configuration?
        '''
        pass

    # knowledge
    PDB_NONSTRU_INFO_LINE_NAME = ['HEADER','TITLE','COMPND','SOURCE','KEYWDS','EXPDTA',
    'NUMMDL','AUTHOR','REVDAT','SPRSDE','JRNL','REMARK','DBREF','SEQADV','SEQRES','HELIX'
    ,'SHEET','CRYST1','ORIGX1','ORIGX2','ORIGX3','SCALE1','SCALE2','SCALE3'] #these keywords are for general information in the PDB database

    # interface
    @classmethod
    def get_structure(cls, path: str, model: int=0,
                           add_solvent_list: List=None, add_ligand_list: List=None, 
                           remove_trash: bool=True, give_idx_map: bool = False) -> Union[Structure, tuple]:
        '''
        Converting a PDB file (as its path) into the Structure()
        Arg:
            path: the file path of the PDB file
            model: The selected model index if there are multiple models in the file
                    (start from 0 default: 0)
                    (assume MODEL appears in order in the file)
            add_solvent_list: (used for categorize residues) additional names for solvent
            add_ligand_list: (used for categorize residues) additional names for ligands
                * solvent list have higher pirority
            remove_trash: if remove trash ligands
            give_idx_map: if return a tuple of (Structure, idx_change_mapper) - which is a dict 
                          with {(old_chain_id, old_residue_id): (new_chain_id, new_residue_id), ... }

        Return:
            Structure()
        '''
        _LOGGER.debug(f'working on {path}')
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
        cls._resolve_missing_chain_id(target_model_df, target_model_ter_df) # add missing chain id in place
        cls._resolve_alt_loc(target_model_df) # resolve alt loc record in place by delele redundant df rows
        #end region (Add here to address more problem related to the PDB file format here)
        
        # target_model_df below should be 'standard' --> start building
        # mapper is for indicating superior information of the current level: e.g.: indicate
        # which residue and chain the atom belongs to in atom_mapper. Another workaround in
        # old enzyhtp is build iteratively chain frist and residues are build in the chain
        # builder and so on for atom. But current methods is used for better readability
        atom_mapper: Dict[tuple, Atom] = cls._build_atoms(target_model_df)
        res_mapper, idx_change_mapper = cls._build_residues(atom_mapper, add_solvent_list, add_ligand_list)
        chain_list: Dict[str, Chain] = cls._build_chains(res_mapper, remove_trash)

        if give_idx_map:
            return Structure(chain_list), idx_change_mapper
        return Structure(chain_list)

    @classmethod
    def get_file_str(cls, stru: Structure) -> str:
        '''
        Convert Structure() into PDB file string
        '''
        pass

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
                _LOGGER.error(
                    f"The PDB '{pdb_path}' contains non-ASCII text and is invalid in line {idx}: '{ll}'. Exiting..."
                )
                sys.exit(1)
   
    @staticmethod
    def _get_target_model(df: pd.DataFrame, model: int) -> Union[pd.DataFrame, None]:
        '''
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
        '''
        if 'MODEL' in df['OTHERS']['record_name'].values:
            # san check for start/end indicator
            mdl_ind_sorted_w_lines = df['OTHERS'][(df['OTHERS'].record_name == 'MODEL') | (df['OTHERS'].record_name == 'ENDMDL')].sort_values(by='line_idx', inplace=False)
            # check if they appear alternately
            first = 1
            for ind in mdl_ind_sorted_w_lines.values:
                if first:
                    last_ind = ind
                    first = 0
                    continue
                if ind[0] == last_ind[0]:
                    _LOGGER.error(f'Unclosed MODEL section is detected in pdb file. line: {last_ind[-1]}')
                    sys.exit(1)
                last_ind = ind
            # get target model unit range
            mdl_start_lines = list(df['OTHERS'][df['OTHERS'].record_name == 'MODEL']['line_idx'])
            mdl_end_lines = list(df['OTHERS'][df['OTHERS'].record_name == 'ENDMDL']['line_idx'])
            mdl_range = list(zip(mdl_start_lines, mdl_end_lines))
            target_mdl_range = mdl_range[model]
            mdl_query_pattern = f'line_idx > {target_mdl_range[0]} & line_idx < {target_mdl_range[1]}'
            # get model dataframe section as a copy
            target_mdl_df = pd.concat((df['ATOM'], df['HETATM'])).query(mdl_query_pattern).copy()
            target_mdl_ter_df = df['OTHERS'][df['OTHERS'].record_name == 'TER'].query(mdl_query_pattern).copy()
        else:
            # get all dataframe as a copy if there's no MODEL record
            target_mdl_df = pd.concat((df['ATOM'], df['HETATM']), ignore_index=True).copy() # insure unique df id
            target_mdl_ter_df = df['OTHERS'][df['OTHERS'].record_name == 'TER'].copy()

        return target_mdl_df, target_mdl_ter_df

    @classmethod
    def _resolve_missing_chain_id(cls, df: pd.DataFrame, ter_df: pd.DataFrame) -> None:
        '''
        Function takes the dataframes of chains and ensures consistent naming
        of chains with no blanks.

        add missing chain id into the Atomic 'df' when there is none
        according to the TER record provide by 'ter_df'
        * ter df should share the same line index references as df (from same file)
        [edit df in place]
        '''
        # no problem case
        chain_ids = set(list(map(lambda c_id: c_id.strip(), df['chain_id'])))
        if '' not in chain_ids:
            return df

        # problem case
        chain_ids.remove('')
        legal_ids = cls._get_legal_pdb_chain_ids(chain_ids)
        ## split the df into chain dfs
        ter_line_ids = list(ter_df['line_idx'])
        chain_dfs = split_df_base_on_column_value(df, 'line_idx', ter_line_ids)
        ## get a list of (loc, value) for editing in the original df
        result_loc_map = []
        for chain in chain_dfs:
            # check for redundant chains
            if len(chain) == 0:
                _LOGGER.debug('Found empty chain. ')
                continue
            atom_missing_c_id = chain[chain.chain_id.str.strip() == '']
            # check if not missing any chain_id
            if len(atom_missing_c_id) == 0:
                    continue
            # check existing chain ids in the chain
            existing_chain_ids = set(list(map(lambda c_id: c_id.strip(), chain['chain_id'])))
            existing_chain_ids.discard('')
            if len(existing_chain_ids) > 1:
                _LOGGER.error('Found multiple chain id in 1 chain. Check your PDB.')
                sys.exit(1)
            # case: all atoms are missing chain id
            if len(existing_chain_ids) == 0:
                result_loc_map.extend(list(zip(atom_missing_c_id.index, [legal_ids.pop()]*len(atom_missing_c_id))))
            # case: only some atoms are missing chain id
            else:
                current_chain_id = list(existing_chain_ids)[0]
                result_loc_map.extend(list(zip(atom_missing_c_id.index, [current_chain_id]*len(atom_missing_c_id))))
        # add missing chain id
        batch_edit_df_loc_value(df, result_loc_map, 'chain_id')

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
    def _resolve_alt_loc(df: pd.DataFrame, keep: str = 'first') -> None:
        '''
        Resolves atoms with the alt_loc records
        Optional argument "keep" specifies the resolution method. Default
        is "first" which forces keeping atoms in the frist residue (earlist 
        in the file) out of multiple ones. Otherwise, keeps the specific 
        alt_loc identifier (e.g.: B), assuming it exists.

        Only one record out of multiple alt_loc is allowed. Delete rest 
        df lines in place. 
        TODO support residue specific keep strageties
        '''
        alt_loc_atoms_df = df[df['alt_loc'].str.strip() != '']
        # san check
        if len(alt_loc_atoms_df) == 0:
            _LOGGER.debug('No alt_loc to resolve.')
            return
        # solve
        # get a list of 'loc' for deleting in the original df
        delete_loc_list = []
        # treat in residues
        alt_loc_residues = alt_loc_atoms_df.groupby(['residue_number','chain_id'], sort=False)
        for r_id_c_id, res_df in alt_loc_residues:
            alt_res_dfs = res_df.groupby('alt_loc', sort=False)
            if len(alt_res_dfs) == 1:
                _LOGGER.debug(f'Only 1 alt_loc id found in residue {r_id_c_id[1], r_id_c_id[0]}. No need to resolve')
                continue
            _LOGGER.debug(f'Dealing with alt_loc residue: {r_id_c_id[1], r_id_c_id[0]}')
            if keep == 'first':
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
        _LOGGER.debug(f'deleting df row: {delete_loc_list}')
        df.drop(index = delete_loc_list, inplace=True)

    @staticmethod
    def _build_atoms(df: pd.DataFrame) -> Dict[str, Atom]:
        '''
        build atom objects from a 'standard' dataframe
        * HETATM should be in seperate chains! HETATM was meant for atoms in the small molecule
          so they cannot be in the same chain as the ATOM atom https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        Return:
            {('chain_id','residue_idx', 'residue_name', 'record_name') : [Atom_obj, ...], ...}
        '''
        atom_mapper = defaultdict(list)
        for i, row in df.iterrows():
            atom_obj = Atom(row)
            residue_key = (row['chain_id'].strip(), row['residue_number'], row['residue_name'].strip(), row['record_name'].strip())
            atom_mapper[residue_key].append(atom_obj)
        return atom_mapper

    @classmethod
    def _build_residues(cls, atom_mapper: Dict[str, Atom], add_solvent_list: List =None, add_ligand_list: List=None) -> Dict[str, Residue]: #TODO add test for this
        '''
        build Residue() objects from atom_mapper
        Arg:
            atom_mapper
            add_solvent_list: (used for categorize residues) additional names for solvent
            add_ligand_list: (used for categorize residues) additional names for ligands
            * solvent list have higher pirority
        Return:
            {'chain_id' : [Residue, ...]}
        '''
        build_mapper = defaultdict(lambda: defaultdict(list))
        result_mapper = defaultdict(list)
        idx_change_mapper = {}
        for res_key, atoms in atom_mapper.items():
            res_obj = Residue(int(res_key[1]), res_key[2] ,atoms)
            build_mapper[res_key[0]][res_key[3]].append(res_obj)
        # give HETATM residue individual chain
        legal_chain_ids = cls._get_legal_pdb_chain_ids(build_mapper.keys())
        for chain_id, record_residues in build_mapper.items():
            if len(record_residues) > 1:
                new_chain_id = legal_chain_ids.pop()
                _LOGGER.debug(f'found HETATM in a ATOM chain: making a new chain {new_chain_id}') #TODO add a recorder for the index change
                result_mapper[new_chain_id] = record_residues['HETATM']
                for i in record_residues['HETATM']:
                    idx_change_mapper[(chain_id, i.idx)] = (new_chain_id, i.idx)
            # in all cases keep the ATOM chain
            result_mapper[chain_id] = list(record_residues.values())[0]
        # categorize_residue
        cls._categorize_pdb_residue(result_mapper, add_solvent_list, add_ligand_list)

        return result_mapper, idx_change_mapper

    @staticmethod
    def _categorize_pdb_residue(residue_mapper: Dict[str, List[Residue]], add_solvent_list: list = None, add_ligand_list: list = None) -> Union[Residue, Ligand, Solvent, MetalUnit]: #TODO add test for this
        '''
        Categorize Residue base on it 3-letter name and chain info in PDB format
        Takes a mapper of {chain_id, [Residue, ...]} and converts
        them into its specialized Residue() inherited class.
        '''
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
                            _LOGGER.error(f'a metal or solvent residue name is found in an peptide chain {chain_id}: {residue.idx} {residue.name}')
                            sys.exit(1)
                        # TODO maybe a class for non-canonical aa in the future
                        _LOGGER.debug(f'found noncanonical {chain_id} {residue.idx}')
                        residue.rtype = chem.ResidueType.NONCANONICAL
                continue
            # non-peptide chain
            for i, residue in enumerate(residues):
                # if metal
                if residue.name in chem.METAL_MAPPER:
                    _LOGGER.debug(f'found metal {chain_id} {residue.idx}')
                    residue_mapper[chain_id][i] = residue_to_metal(residue)
                    continue
                # if solvent
                if residue.name in chem.RD_SOLVENT_LIST + add_solvent_list:
                    _LOGGER.debug(f'found solvent {chain_id} {residue.idx}')
                    residue_mapper[chain_id][i] = residue_to_solvent(residue)
                    continue
                if residue.name in (set(chem.RD_NON_LIGAND_LIST) - set(add_ligand_list)):
                    _LOGGER.debug(f'found trash {chain_id} {residue.idx}')
                    residue.rtype = chem.ResidueType.TRASH
                    continue
                # the left cases are ligands
                _LOGGER.debug(f'found ligand {chain_id} {residue.idx}')
                residue_mapper[chain_id][i] = residue_to_ligand(residue)

    @staticmethod
    def _build_chains(res_mapper: Dict[str, Residue], remove_trash: bool) -> Dict[str, Chain]:
        '''
        build Chain() objects from residue mapper
        '''
        chain_list = []
        for chain_id, residues in res_mapper.items():
            ch_obj = Chain(chain_id ,residues) # TODO decouple chain name with pdb chain id? How to solve too many water chain
            if remove_trash:
                ch_obj.remove_trash()
            chain_list.append(ch_obj)
        return chain_list
    #endregion

    @staticmethod
    def _make_pdb_atom():
        pass


# TODO go to core helper
def split_df_base_on_column_value(df: pd.DataFrame, column_name: str, split_values: list, copy: bool=False):
    '''
    split a dataframe base on the value of a column
    ** the line in the split values will not be included **
    Arg:
        df: the target dataframe
        column_name: the reference column's name
        split_value: the value mark for spliting
    '''
    split_values = sorted(split_values)
    frist = 1
    result_dfs = []
    for this_value in split_values:
        if frist:
            frist_df = df[df[column_name] < this_value]
            result_dfs.append(frist_df)
            frist = 0
            last_value = this_value
            continue
        result_df = df[(df[column_name] < this_value) & (df[column_name] > last_value)]
        result_dfs.append(result_df)
        last_value = this_value
    # deal with the last portion
    result_dfs.append(df[df[column_name] > split_values[-1]])

    if copy:
        for i, df in enumerate(result_dfs):
            result_dfs[i] = df.copy()

    return result_dfs

def batch_edit_df_loc_value(df: pd.DataFrame, loc_value_list: List[tuple], column: str):
    '''
    batch edit 'column' of 'df' with the 'loc_value_list'
    '''
    for loc, value in loc_value_list:
        df.loc[loc, column] = value