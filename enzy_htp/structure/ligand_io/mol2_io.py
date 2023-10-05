"""Generation/construction of Ligand() objects from .mol2 files and exporting these objects to this file format. 
Definition of .mol2 format (http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf). All parsing is done within enzy_htp using 
this parser only. The Mol2Parser has no private data and serves as a namespace for .mol2 I/O conversion functions.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-10-05
"""
from typing import List, Dict, Any

import pandas as pd


from enzy_htp.core import _LOGGER
import enzy_htp.core.file_system as fs

from ..atom import Atom
from ..ligand import Ligand

class Mol2Parser():
    """Holds all functionality for .mol2 I/O with respect to the Ligand() class. This parser contains no private data
    and instead serves as a namespace to hold this functionality. All functions are static classmethods. Functionality is intended to 
    be accessed through the get_ligand() and save_ligand() methods. All other methods are meant for implementation only and SHOULD NOT
    be accessed directly by users. 
    """
   

    @classmethod
    def _parse_molecule_section(cls, raw:List[str]) -> Dict:
        """Parses out the @<TRIPOS>MOLECULE section of a .mol2 file and returns the contents in a dictionary. Keys in the dictionary
        can include:
            mol_name, mol_type, charge_type, num_atoms, num_bond, num_subst, num_feat, num_sets

        Args:
            raw: List[str] of the raw contents.

        Returns:
            A dictionary with the above keys. 

        """
        result:Dict = {
            'mol_name':raw[0].strip(),
            'mol_type':raw[2].strip(),
            'charge_type':raw[3].strip()
        }

        for label, raw_value in zip("num_atoms num_bond num_subst num_feat num_sets".split(), raw[1].split()):
            result[label] = int(raw_value)
    
        return result

    @classmethod
    def _parse_atom_section(cls, raw:List[str] ) -> List[Dict]:
        """Parses out the @<TRIPOS>ATOM section of a .mol2 file and returns the contents as a List[Dict] where each Dict
        represents the information in a single line, representing a single atom. Keys in the dictionary can include:
            atom_id atom_name x_coord y_coord z_coord atom_type subst_id subst_name charge status_bit

        Args:
            raw: List[str] of the raw contents.

        Returns:
            A List[Dict] which has one entry for each atom in the .mol2 file.
        """
        result = list()
        labels:List[str]="atom_id atom_name x_coord y_coord z_coord atom_type subst_id subst_name charge status_bit".split()
        types:List[Any] =[    int,      str,  float,  float,  float,      str,     int,       str, float,       str]
        for rr in raw:
            temp = dict()
            for ll, tt, cc in zip(labels, rr.split(), types):
                temp[ll] = cc(tt)
            result.append(temp)                

        return result

    @classmethod
    def _parse_bond_section(cls, raw:List[str] ) -> List[Dict]:
        """Parses out the @<TRIPOS>BOND section of a .mol2 file and returns the contents as a List[Dict] where each Dict
        represents the information in a single line, representing a single bond. Keys in the dictionary can include:
            bond_id origin_atom_id target_atom_id bond_type status_bits

        Args:
            raw: List[str] of the raw contents.

        Returns:
            A List[Dict] which has one entry for each bond in the .mol2 file.
        """
        result = list()
        labels:List[str]="bond_id origin_atom_id target_atom_id bond_type status_bits".split()
        types:List[Any] =[    int,           int,           int,      str,         str]
        for rr in raw:
            temp = dict()
            for ll, tt, cc in zip(labels, rr.split(), types):
                temp[ll] = cc(tt)
            result.append(temp)

        return result 

    @classmethod
    def _parse_substructure_section(cls, raw:List[str]) -> List[Dict]:
        """Parses out the @<TRIPOS>SUBSTRUCTURE section of a .mol2 file and returns the contents as a List[Dict] where each Dict
        represents the information in a single line, representing a single substructure. Keys in the dictionary can include:
            subst_id subst_name root_atom subst_type dict_type chain sub_type inter_bonds status

        Args:
            raw: List[str] of the raw contents.

        Returns:
            A List[Dict] which has one entry for each substructure in the .mol2 file.
        """

        result = list()
        labels:List[str]="subst_id subst_name root_atom subst_type dict_type chain sub_type inter_bonds status".split()
        types:List[Any] =[     int,       str,      int,       str,      str,  str,     str,        int,   str]
        for rr in raw:
            temp = dict()
            for ll, tt, cc in zip(labels, rr.split(), types):
                temp[ll] = cc(tt)
            result.append( temp )

        return result

    @classmethod
    def parse_mol2_file(cls, path:str) -> Dict[str,List[Any]]:
        """Fast and simple parser for .mol2 file format. All information is stored in and returned in a dict() where
        the keys are @<TRIPOS> headers ("ATOM", "BOND", "MOLECULE", etc.). The values for each of the keys may be 
        parsed into a Dict or List[Dict], but may also be the original, raw List[str] from the input file. Function also
        performs basic checks, ensuring that the number of bonds and atoms match the expected values. Also checks that the
        input file exists and is not empty.
        
        Args:
            path: A str() path to a .mol2 file to parse.

        Returns:
            A Dict with (key, value) pairs that mirror the @<TRIPOS> headings in the file.

        """
        if not fs.has_content( path ):
            _LOGGER.error(f"The supplied file {path} does not exist or is empty. Exiting...")
            exit( 1 )
                
        result = dict()
        content:List[str]=fs.lines_from_file( path )
        content = list(filter(lambda ll: not cls.is_comment( ll ), content ))
        for cidx,cc in enumerate(content):
            if cc.startswith('@<TRIPOS>'):
                it:int = cidx+1
                temp = list()
                while it < len(content) and content[it].find('@<TRIPOS>') == -1:
                    temp.append(content[it])
                    it += 1
                result[cc.replace('@<TRIPOS>','')] = temp

        for k,v in result.items():
            if k == 'MOLECULE':
                result[k] = cls._parse_molecule_section(v)
            elif k == 'BOND':
                result[k] = cls._parse_bond_section(v)
            elif k == 'ATOM':
                result[k] = cls._parse_atom_section(v)
            elif k == 'SUBSTRUCTURE':
                result[k] = cls._parse_substructure_section(v)

        
        error:bool = False
        if result['MOLECULE']['num_atoms'] != len(result['ATOM']):
            error = True
            _LOGGER.error(f"The number of atoms {len(result['ATOM'])} in {path} is not consistent with the expected amount number {result['MOLECULE']['num_atoms']}")

        if result['MOLECULE']['num_bond'] != len(result['BOND']):
            error = True
            _LOGGER.error(f"The number of bonds {len(result['BOND'])} in {path} is not consistent with the expected amount number {result['MOLECULE']['num_bond']}")


        if error:
            _LOGGER.error("Exiting...")
            exit( 1 )

        return result

    @classmethod
    def is_comment(cls, line :str ) -> bool:
        """Is the supplied line a comment as defined in the .mol2 format?"""
        tks = line.split('#')
        return not len(tks[0])
        

    @classmethod
    def get_ligand(cls,
                        path :str,
                        residue_name:str=None,
                        residue_idx:int=None,
                        chain_name:str=None
                        ) -> Ligand:
        """Creates a Ligand() object from a .mol2 file. Assumes that the input file has exactly one structure and has
        undefined behavior if this is not the case. Primarily uses ligand, atom, and bond information to define new Ligand() object.

        Args:
            path: The filename to load the Ligand() from as a str().
            residue_name: The 3-letter PDB style name of the residue as a str().
            residue_idx: The int() index of the Ligand().
            chain_name: The str() chain name of the Ligand().

        Returns:
            The newly created Ligand() object.
        """

        file_info:Dict=cls.parse_mol2_file(path)

        subst_name:str = file_info['SUBSTRUCTURE'][0]['subst_name']

        if residue_name is None:
            residue_name = subst_name[0:3]

        if residue_idx is None:
            residue_idx = int(subst_name[3:])

        atoms = list()
        for aa in file_info['ATOM']:
            labels:List[str]="x_coord y_coord z_coord atom_name charge atom_type".split()
            temp = dict()
            for ll in labels:
                temp[ll] = aa[ll]
            atoms.append( Atom(temp) )

        return Ligand(residue_idx=residue_idx, residue_name=residue_name, atoms=atoms, bonds=file_info['BOND'])            

    @classmethod
    def save_ligand(cls, outfile:str, ligand:Ligand) -> str:
        """Given the name of an output .mol2 file and a Ligand() object, saves that molecule to that file. Even if the Ligand()
        is part of a chain, only the Ligand() is saved. Will not work with Residue-derived objects that do not have explicit bonds.

        Args:
            outfile: The name of the .mol2 file to save the Ligand() as a str().
            ligand: The Ligand() object to save.

        Returns:
            The name of the newly saved .mol2 file as a str().
        """
        content:List[str] = [
            "@<TRIPOS>MOLECULE",
            ligand.name if ligand.name else "UNK",
            f"{len(ligand.atoms)} {len(ligand.bonds)} 1",
            "SMALL",
            "USER_CHARGES",
            "@<TRIPOS>ATOM"
        ]
        for aidx,atom in enumerate(ligand.atoms):
            content.append(f"{aidx+1}\t{atom.name: >4}\t{atom.coord[0]:.3f}\t{atom.coord[1]:.3f}\t{atom.coord[2]:.3f}\t{atom._atom_type}\t{1}\t{ligand.name}{ligand.idx}\t{atom.charge:.3f}")

        content.append("@<TRIPOS>BOND")

        for bond in ligand.bonds:
            content.append(f"{bond['bond_id']} {bond['origin_atom_id']} {bond['target_atom_id']} {bond['bond_type']}")

        chain = "****"
        if ligand.chain is not None:
            chain = ligand.chain.name

        content.extend([
            "@<TRIPOS>SUBSTRUCTURE",
           f"1\t{ligand.name}{ligand.idx}\t1\tGROUP\t1 {chain}\t{ligand.name}"
        ])
        fs.write_lines( outfile, content )

        return outfile
            
