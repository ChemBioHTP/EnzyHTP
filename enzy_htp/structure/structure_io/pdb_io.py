"""Generation/construction of Structure objects from PDB files and exporting it vice versa
Definition of PDB file format (http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)

Author: shaoqz, <shaoqz@icloud.com>
Date: 2022-08-01
"""
import os
import string
import sys
from typing import Dict, List, Tuple, Union
from biopandas.pdb import PandasPdb
import pandas as pd

from ._interface import StructureParserInterface
from ..atom import Atom
from ..metal_atom import MetalAtom, residue_to_metal
from ..residue import Residue
from ..solvent import Solvent, residue_to_solvent
from ..ligand import Ligand, residue_to_ligand
from ..chain import Chain
from ..structure import Structure
from enzy_htp.core import _LOGGER
import enzy_htp.core.file_system as fs

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
    def get_structure(cls, path: str, model: int=0) -> Structure:
        '''
        Converting a PDB file (as its path) into the Structure()
        Arg:
            path: the file path of the PDB file
            model: The selected model index if there are multiple models in the file
                    (start from 0 default: 0)
                    (assume MODEL appears in order in the file)
        '''
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
        cls._resolve_alt_loc(target_model_df) # resolve alt loc record in place
        #end region (Add here to address more problem related to the PDB file format here)
        # target_model_df below should be 'standard'

        # # mapper is for indicating superior information of the current level: e.g.: indicate
        # # which residue and chain the atom belongs to in atom_mapper. Another workaround in
        # # old enzyhtp is build iteratively chain frist and residues are build in the chain
        # # builder and so on for atom. But current methods is used for better readability
        # atom_mapper: Dict[str, Atom] = cls._build_atoms(target_model_df) 
        # res_mapper: Dict[str, Residue] = cls._build_residues(atom_mapper)
        # chain_list: Dict[str, Chain] = cls._build_chains(res_mapper)

        # return Structure(chain_list)

    @classmethod
    def get_file_str(cls, stru: Structure) -> str:
        '''
        Convert Structure() into PDB file string
        '''
        pass

    # region === PDB -> Stru ===
    @staticmethod
    def _check_valid_pdb(pdb_path: str) -> None:
        """Helper function that ensures the supplied pdb_path contain a valid pdb file.
        Private to structure.structure_io.pdb_io.py. Should NOT be called externally."""
        # pylint: disable=logging-fstring-interpolation, consider-using-sys-exit
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
            target_mdl_df = pd.concat((df['ATOM'], df['HETATM'])).copy()
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
        ## get a list of (loc, value) part for editing in the original df
        result_loc_map = [] 
        for chain in chain_dfs:
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
    def _get_legal_pdb_chain_ids(taken_ids: List[str]) -> List[str]:
        """
        Small helper method that determines the legal PDB chain names given a list of taken ids.
        Uses all of the available 26 capitalized letters and returns in reverse order. Returns an 
        empty list when all are occupied.
        TODO check how PDB defines 26+ chains (maybe AA AB etc.)
        Return:
            the legal chain name list in a reversed order (for doing .pop())
        """
        result = list(string.ascii_uppercase)
        result = list(filter(lambda s: s not in taken_ids, result))
        return list(reversed(result))

    @staticmethod
    def _resolve_alt_loc():
        pass

    @staticmethod
    def _build_atom(df: pd.DataFrame):
        '''
        create atom objects from a dataframe
        '''
        pass

    @staticmethod
    def _build_residues():
        pass

    @staticmethod
    def _build_chains():
        pass
    # endregion

class PDBAtom():
    '''
    class for accessing atom records in the PDB file
    contructed with a pandas.Series
    '''
    pass

# TODO go to core helper
def split_df_base_on_column_value(df: pd.DataFrame, column_name: str, split_values: list, copy: bool=False):
    '''
    split a dataframe base on the value of a column
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