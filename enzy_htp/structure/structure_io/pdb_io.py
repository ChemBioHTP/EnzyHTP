"""Generation/construction of Structure objects from PDB files and exporting it vice versa
Definition of PDB file format (http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)

Author: shaoqz, <shaoqz@icloud.com>
Date: 2022-08-01
"""
import os
from typing import Dict, Union
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
        input_pdb.read_pdb(path) # this add a dataframe to the obj
        # build Structure object
        
        # deal with multiple model
        target_model_df, target_model_ter_df = cls._get_target_model(input_pdb.df, model=model)
        
        用对应的model建造atom 然后用ter df补全链标
        # atom_mapper: Dict[str, Atom] = build_atom(
        #     pass
        # )
        # res_mapper: Dict[str, Residue] = build_residues(
        #     pd.concat((input_pdb.df["ATOM"], input_pdb.df["HETATM"])), keep #@shaoqz: so it works with multichains w/o chain id? I saw that four chains in test file
        # )
        # chain_mapper: Dict[str, Chain] = build_chains(res_mapper)
        # stru_obj = Structure(list(chain_mapper.values()))
        # return stru_obj


    @classmethod
    def get_file_str(cls, stru: Structure) -> str:
        '''
        Convert Structure() into PDB file string
        '''
        pass

    #region === PDB -> Stru ===
    @staticmethod
    def _check_valid_pdb(pdb_path: str) -> None:
        """Helper function that ensures the supplied pdb_path contain a valid pdb file.
        Private to structure.structure_io.pdb_io.py. Should NOT be called externally."""
        # pylint: disable=logging-fstring-interpolation, consider-using-sys-exit
        # check for right extension
        ext: str = fs.get_file_ext(pdb_path)
        if ext.lower() != ".pdb":
            _LOGGER.error(f"Supplied file '{pdb_path}' is NOT a PDB file. Exiting...")
            exit(1)
        #check for file existance
        if not os.path.exists(pdb_path):
            _LOGGER.error(f"Supplied file '{pdb_path}' does NOT exist. Exiting...")
            exit(1)
        # check for character encoding
        for idx, ll in enumerate(fs.lines_from_file(pdb_path)):
            if not ll.isascii():
                _LOGGER.error(
                    f"The PDB '{pdb_path}' contains non-ASCII text and is invalid in line {idx}: '{ll}'. Exiting..."
                )
                exit(1)
    
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
        '''
        if not 'MODEL' in df['OTHERS']['record_name'].values:
            return

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
                exit(1)
            last_ind = ind
        # get target model unit range
        mdl_start_lines = list(df['OTHERS'][df['OTHERS'].record_name == 'MODEL']['line_idx'])
        mdl_end_lines = list(df['OTHERS'][df['OTHERS'].record_name == 'ENDMDL']['line_idx'])
        mdl_range = list(zip(mdl_start_lines, mdl_end_lines))
        target_mdl_range = mdl_range[model]
        # get model dataframe section
        target_mdl_df = pd.concat((df['ATOM'], df['HETATM'])).query(f'line_idx > {target_mdl_range[0]} & line_idx < {target_mdl_range[1]}')
        target_mdl_ter_df = df['OTHERS'][df['OTHERS'].record_name == 'TER'].query(f'line_idx > {target_mdl_range[0]} & line_idx < {target_mdl_range[1]}')
        return target_mdl_df, target_mdl_ter_df

    # def legal_chain_names(mapper: Dict[str, List[Residue]]) -> List[str]:
    #     """Small helper method that determines the legal chain names given a Residue() mapper.
    #     Uses all of the available 26 capitalized letters and returns in reverse order. Returns an empty list when all are occupied."""
    #     result = list(string.ascii_uppercase)
    #     taken = set(list(mapper.keys()))
    #     result = list(filter(lambda s: s not in taken, result))
    #     return list(reversed(result))


    # def name_chains(mapper: Dict[str, List[Residue]]) -> None:
    #     """Function takes a defaultdict(list) of Residues and ensures consistent naming of chains with no blanks."""
    #     key_names = set(list(map(lambda kk: kk.strip(), mapper.keys())))
    #     if "" not in key_names:
    #         return mapper
    #     unnamed = list(mapper[""])
    #     del mapper[""]

    #     names = legal_chain_names(mapper)
    #     unnamed = sorted(unnamed, key=lambda r: r.min_line())
    #     new_chain: List[Residue] = []

    #     for res in unnamed:
    #         if not new_chain:
    #             new_chain.append(res)
    #         elif new_chain[-1].neighbors(res):
    #             new_chain.append(res)
    #         else:
    #             mapper[names.pop()] = deepcopy(new_chain)
    #             new_chain = [res]
    #     if new_chain:
    #         mapper[names.pop()] = deepcopy(new_chain)
    #     return mapper


    # def categorize_residue(residue: Residue) -> Union[Residue, Ligand, Solvent, MetalAtom]:
    #     """Method that takes a default Residue() and converts it into its specialized Residue() inherited class."""
    #     # TODO(CJ): I need to add in stuff here for the solvent types.
    #     if residue.is_canonical():
    #         residue.set_rtype(chem.ResidueType.CANONICAL)
    #         return residue

    #     if (
    #         residue.name in chem.METAL_MAPPER
    #     ):  # TODO(CJ): implement more OOP method for this
    #         return residue_to_metal(residue)

    #     if residue.is_rd_solvent():  # TODO(CJ): make sure this logic is 100% right
    #         return residue_to_solvent(residue)

    #     return residue_to_ligand(residue)


    # def build_residues(
    #     df: pd.DataFrame, keep: str = "first"
    # ) -> Dict[str, Union[Residue, Ligand, Solvent, MetalAtom]]:
    #     """Helper method that builds Residue() or derived Residue() objects from a dataframe generated by BioPandas.
    #     Returns as a dict() with (key, value) pairs of (residue_key, Union[Residue,Ligand,Solvent,MetalAtom]).""" #@shaoqz: should we also explain the keep here?
    #     mapper = defaultdict(list)
    #     for i, row in df.iterrows():
    #         aa = Atom(**row)
    #         mapper[aa.residue_key()].append(aa)
    #     # for k,v in mapper.items():
    #     # print(k,len(v))
    #     result: Dict[str, Residue] = dict()
    #     for res_key, atoms in mapper.items():
    #         result[res_key] = Residue(
    #             residue_key=res_key, atoms=sorted(atoms, key=lambda a: a.atom_number)
    #         )
    #     for (res_key, res) in result.items():
    #         result[res_key] = categorize_residue(res)
    #         if keep != "all":
    #             result[res_key].resolve_alt_loc(keep)
    #             result[res_key].remove_alt_loc()
    #             # result[res_key].remove_occupancy() #@shaoqz: what happened here?
    #     return result


    # def build_chains(mapper: Dict[str, Residue]) -> Dict[str, Chain]:
    #     """Helper method that builds the Chain() objects from a dict() of (residue_key, Residue()) pairs generated by build_residues()."""
    #     chain_mapper = defaultdict(list)
    #     for res in mapper.values():
    #         chain_mapper[res.chain()].append(res)
    #     chain_mapper = name_chains(chain_mapper)
    #     result: Dict[str, Chain] = dict()
    #     # ok this is where we handle missing chain ids
    #     for chain_name, residues in chain_mapper.items():
    #         result[chain_name] = Chain(chain_name, sorted(residues, key=lambda r: r.num()))
    #     return result

    # def ligand_from_pdb(fname: str, net_charge: float = None) -> Ligand:
    #     """Creates a Ligand() object from a supplied .pdb file. Checks that the input file both exists and is ASCII format."""
    #     def get_charge_mapper(fname:str) -> dict:
    #         """Helper method that gets charges"""
    #         lines = fs.lines_from_file(fname) 
    #         result = dict()
    #         for ll in lines:
    #             raw_charge = ll[78:80].strip()
    #             if not raw_charge:
    #                 continue
    #             mult = 1 
    #             temp = ''
    #             for ch in raw_charge:
    #                 if ch != '-':
    #                     temp += ch
    #                 else:
    #                     mult = -1
    #             result[(ll[21].strip(), int(ll[6:11]))] = mult*int( temp )
    #         return result 

    #     check_valid_pdb(fname)
    #     warnings.filterwarnings("ignore")
    #     # adapt general input // converge to a list of PDB_line (resi_lines)
    #     parser = PandasPdb()
    #     parser.read_pdb(fname)
    #     temp_df = pd.concat((parser.df["ATOM"], parser.df["HETATM"]))
    #     charge_mapper = get_charge_mapper( fname )
    #     #TODO(CJ): put this into its own function
    #     for (cname, a_id), charge in charge_mapper.items():
    #         #print(cname, a_id, charge)
    #         mask = ((temp_df.chain_id==cname)&(temp_df.atom_number==a_id)).to_numpy()
    #         idx = np.where(mask)[0][0]
    #         temp_df.at[idx, 'charge'] = charge
    #     atoms = list(
    #         map(
    #             lambda pr: Atom(**pr[1]),
    #             temp_df.iterrows(),
    #         )
    #     )
    #     residue_name = set(list(map(lambda a: a.residue_name, atoms)))
    #     residue_number = set(list(map(lambda a: a.residue_number, atoms)))
    #     chain_id = set(list(map(lambda a: a.chain_id, atoms)))
    #     assert len(residue_name) == 1
    #     assert len(residue_number) == 1
    #     assert len(chain_id) == 1
    #     # TODO(CJ): should I set the chain for this if it is blank?
    #     key = f"{list(chain_id)[0]}.{list(residue_name)[0]}.{list(residue_number)[0]}"
    #     result = Residue(key, atoms)
    #     result = residue_to_ligand(result)
    #     result.net_charge = net_charge  # TODO(CJ): make net_charge a getter/setter
    #     return result

    # def get_ligand_name(fname: str) -> str:
    #     # TODO(CJ): add the documentation herej
    #     # TODO(CJ): make this more efficient
    #     # TODO(CJ): add testing for this
    #     ligand: Ligand = ligand_from_pdb(fname)
    #     return deepcopy(ligand.get_name())
        
    #endregion
