"""Generation/construction of Structure() objects from .prmtop files and exporting these objects to this file format. 
Definition of .prmtop format (http://ambermd.org/FileFormats.php#topology). All parsing is done within enzy_htp using 
this parser only. The PrmtopParser has no private data and serves as a namespace for .prmtop I/O conversion functions.

Author: Qianzhen Shao <shaoqz@icloud.com>
Date: 2023-10-28
"""
from collections import defaultdict
from typing import Dict, List
import re
from itertools import chain
import periodictable

from ._interface import StructureParserInterface
from ..structure import Structure, Atom, Residue, Chain
from .pdb_io import PDBParser

from enzy_htp.chemical.solvent import RD_SOLVENT_LIST
import enzy_htp.core.file_system as fs
import enzy_htp.core.fortran_helper as fh
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import FileFormatError

class PrmtopParser(StructureParserInterface):
    """the parser for Gaussian prmtop files"""
    def __init__(self, ncaa_chrgspin_mapper: Dict = None) -> None:  # pylint: disable=super-init-not-called
        if ncaa_chrgspin_mapper is None:
            ncaa_chrgspin_mapper = {}
        self.ncaa_chrgspin_mapper = ncaa_chrgspin_mapper

    def get_structure(
            self,
            path: str,
            add_solvent_list: List = None,
            add_ligand_list: List = None,
        ) -> Structure:
        """Converting a .prmtop file (as its path) into the Structure()
        also assign charges as prmtop contains.
        NOTE that this structure dont have actual coordinates/geometry
        Arg:
            path:
                the file path of the prmtop file
            add_solvent_list:
                (used for categorize residues) additional 
                solvents are recognized by matching names 
                non-polypeptide chains:
                    chem.RD_SOLVENT_LIST + add_solvent_lis
            add_ligand_list:
                (used for categorize residues) additional 
                ligands are recognized by not matching nam
                in non-polypeptide chains:
                    1. chem.RD_SOLVENT_LIST + add_solvent_
                    2. chem.RD_NON_LIGAND_LIST - add_ligan
                    * solvent list have higher pirority
        Return:
            Structure()"""
        _LOGGER.debug(f"Building Structure() from prmtop ({path}) ")
        if add_solvent_list is None:
            add_solvent_list = []
        data = type(self)._parse_prmtop_file(path)

        # san check
        essential = ["ATOM_NUMBER", "RESIDUE_CHAINID", "RESIDUE_NUMBER"]
        for ess in essential:
            if ess not in data:
                _LOGGER.error(f"{ess} not in prmtop. ({path}) You need a prmtop after `add_pdb`.")
                raise ValueError
        
        # build atoms (contain charge)
        atoms = []
        atom_data_list = zip(data["ATOM_NAME"], data["CHARGE"], data["ATOMIC_NUMBER"], data["ATOM_NUMBER"])
        for name, charge, ele_num, idx in atom_data_list:
            charge = charge / 18.2223 # unit: e
            atom = Atom(
                name = name,
                coord = (None, None, None),
                idx = idx,
                element = str(periodictable.elements[ele_num]),
                charge = charge,
            )
            atoms.append(atom)

        # build residues
        residue_mapper = defaultdict(list)
        next_pointers = data["RESIDUE_POINTER"][1:] + [None]
        residue_data_list = zip(data["RESIDUE_LABEL"], data["RESIDUE_POINTER"], next_pointers, data["RESIDUE_CHAINID"], data["RESIDUE_NUMBER"])
        taken_ch_ids = set(data["RESIDUE_CHAINID"])
        legal_ch_ids = PDBParser._get_legal_pdb_chain_ids(taken_ch_ids)
        legal_res_idx = type(self)._get_legal_residue_idx(data["RESIDUE_NUMBER"])
        solvent_list = ["Na+", "Cl-"] + RD_SOLVENT_LIST + add_solvent_list
        for name, pointer, next_pointer, chain_id, idx in residue_data_list:
            # resolve chain id
            if chain_id is None:
                if name in solvent_list:
                    chain_id = legal_ch_ids.pop()
                    _LOGGER.debug("Found solvent with out chain id in prmtop. "
                                  f"Assigning new chain: {chain_id}.")
                else:
                    _LOGGER.error("chain id lost in prmtop! "
                                  f"residue {name} from {path} dont have a chain id.")
                    raise ValueError
            # resolve pointer
            res_atoms = type(self)._resolve_residue_pointer(pointer, next_pointer, atoms)
            # resolve missing residue index
            if idx == 0:
                if name in solvent_list:
                    idx = legal_res_idx.pop()
                else:
                    _LOGGER.warning(f"found non-solvent residue index 0 from prmtop. make sure it is correct. ({path})")
            
            residue = Residue(
                residue_name=name,
                residue_idx=idx,
                atoms = res_atoms,
            )
            residue_mapper[chain_id].append(residue)  # here it is {"chain_id": Residue()}

        # categorize_residue
        PDBParser._categorize_pdb_residue(residue_mapper, add_solvent_list, add_ligand_list)

        # build chains
        chains = PDBParser._build_chains(residue_mapper, False)

        # build structure
        result = Structure(chains)
        result.assign_ncaa_chargespin(self.ncaa_chrgspin_mapper)

        return result

    @classmethod
    def _resolve_residue_pointer(cls, pointer: int, next_pointer: int, atoms: List[Atom]) -> List[Atom]:
        """resolve residue pointers from prmtop
        return a list of Atom()s of indicated by the pointer"""
        pointer = pointer - 1
        if next_pointer is None:
            return atoms[pointer:]
        else:
            next_pointer = next_pointer - 1
            return atoms[pointer:next_pointer]
 
    @classmethod
    def _get_legal_residue_idx(cls, taken_ids: List) -> List:
        """get legal residue indexes"""
        max_id = max(taken_ids)
        result = list(range(max_id+1, max_id+10000))
        return list(reversed(result))
    
    @classmethod
    def get_file_str(cls, stru: Structure) -> str:
        """convert a Structure() to .prmtop file content."""
        pass # TODO: add when need

    @classmethod
    def save_structure(cls, out_path: str, stru: Structure) -> str:
        """save a Structure() to .prmtop file.
        return the out_path"""
        pass # TODO: add when need

    @classmethod
    def _parse_prmtop_file(cls, path: str) -> Dict:
        """parse prmtop file to a data dictionary.

        the return data dictionary strictly follows keywords defined
        in https://ambermd.org/prmtop.pdf.
        Example:
            {
            "VERSION_STAMP" : "V0001.000",
            "DATE" : "01/03/24 16:27:03",
            "TITLE" : "default_name",
            ...
            }"""
        if not fs.has_content( path ):
            _LOGGER.error(f"The supplied file {path} does not exist or is empty. Exiting...")
            raise ValueError

        result = {}
        content: str = fs.content_from_file( path )
        sections = re.split(r"^%FLAG", content, flags=re.MULTILINE)

        for sec in sections:
            if sec.find("%VERSION") != -1:
                result.update(cls._parse_version(sec))
                continue
            sec_name = sec.splitlines()[0].strip()
            result.update(cls.section_parser_mapper()[sec_name](sec))

        return result

    @classmethod
    def section_parser_mapper(cls):
        """the mapper that map section type/name to the parser function"""
        return {
            "TITLE" : cls._parse_title,
            "POINTERS" : cls._parse_pointers,
            "ATOM_NAME" : cls._parse_atom_name,
            "CHARGE" : cls._parse_charge,
            "ATOMIC_NUMBER" : cls._parse_atomic_number,
            "MASS" : cls._parse_mass,
            "ATOM_TYPE_INDEX" : cls._parse_atom_type_index,
            "NUMBER_EXCLUDED_ATOMS" : cls._parse_number_excluded_atoms,
            "NONBONDED_PARM_INDEX" : cls._parse_nonbonded_parm_index,
            "RESIDUE_LABEL" : cls._parse_residue_label,
            "RESIDUE_POINTER" : cls._parse_residue_pointer,
            "BOND_FORCE_CONSTANT" : cls._parse_bond_force_constant,
            "BOND_EQUIL_VALUE" : cls._parse_bond_equil_value,
            "ANGLE_FORCE_CONSTANT" : cls._parse_angle_force_constant,
            "ANGLE_EQUIL_VALUE" : cls._parse_angle_equil_value,
            "DIHEDRAL_FORCE_CONSTANT" : cls._parse_dihedral_force_constant,
            "DIHEDRAL_PERIODICITY" : cls._parse_dihedral_periodicity,
            "DIHEDRAL_PHASE" : cls._parse_dihedral_phase,
            "SCEE_SCALE_FACTOR" : cls._parse_scee_scale_factor,
            "SCNB_SCALE_FACTOR" : cls._parse_scnb_scale_factor,
            "SOLTY" : cls._parse_solty,
            "LENNARD_JONES_ACOEF" : cls._parse_lennard_jones_acoef,
            "LENNARD_JONES_BCOEF" : cls._parse_lennard_jones_bcoef,
            "BONDS_INC_HYDROGEN" : cls._parse_bonds_inc_hydrogen,
            "BONDS_WITHOUT_HYDROGEN" : cls._parse_bonds_without_hydrogen,
            "ANGLES_INC_HYDROGEN" : cls._parse_angles_inc_hydrogen,
            "ANGLES_WITHOUT_HYDROGEN" : cls._parse_angles_without_hydrogen,
            "DIHEDRALS_INC_HYDROGEN" : cls._parse_dihedrals_inc_hydrogen,
            "DIHEDRALS_WITHOUT_HYDROGEN" : cls._parse_dihedrals_without_hydrogen,
            "EXCLUDED_ATOMS_LIST" : cls._parse_excluded_atoms_list,
            "HBOND_ACOEF" : cls._parse_hbond_acoef,
            "HBOND_BCOEF" : cls._parse_hbond_bcoef,
            "HBCUT" : cls._parse_hbcut,
            "AMBER_ATOM_TYPE" : cls._parse_amber_atom_type,
            "TREE_CHAIN_CLASSIFICATION" : cls._parse_tree_chain_classification,
            "JOIN_ARRAY" : cls._parse_join_array,
            "IROTAT" : cls._parse_irotat,
            "RADIUS_SET" : cls._parse_radius_set,
            "RADII" : cls._parse_radii,
            "SCREEN" : cls._parse_screen,
            "SOLVENT_POINTERS" : cls._parse_solvent_pointers,
            "ATOMS_PER_MOLECULE" : cls._parse_atoms_per_molecule,
            "BOX_DIMENSIONS" : cls._parse_box_dimensions,
            "CAP_INFO" : cls._parse_cap_info,
            "CAP_INFO2" : cls._parse_cap_info2,
            "POLARIZABILITY" : cls._parse_polarizability,
            "IPOL" : cls._parse_ipol,
            "ATOM_NUMBER" : cls._parse_atom_number,
            "ATOM_BFACTOR" : cls._parse_atom_bfactor,
            "ATOM_OCCUPANCY" : cls._parse_atom_occupancy,
            "ATOM_ELEMENT" : cls._parse_atom_element,
            "RESIDUE_CHAINID" : cls._parse_residue_chainid,
            "RESIDUE_NUMBER" : cls._parse_residue_number,
        }

    @classmethod
    def _parse_version(cls, content: str) -> Dict:
        """parse the 'version' section from prmtop file format
        see test for example of the return."""
        pattern = r"%VERSION  VERSION_STAMP = (.+)DATE = (.+)"
        version_stamp, data = re.match(pattern, content.strip()).groups()

        return {
            "VERSION_STAMP" : version_stamp.strip(),
            "DATE" : data.strip(),            
        }

    @classmethod
    def _parse_general(cls, content: str) -> Dict:
        """parse general data section from prmtop file format"""
        (key, fmt, body) = content.split('\n', 2) # TODO this not work when multiple FORMAT is in the file but it is a rare case when IFPERT == 1
        key = key.strip()
        fmt = fmt.strip().lstrip("%")
        body = re.sub('\n', '', body)
        data = fh.parse_data(fmt, body, reduce=False)
        data = list(chain.from_iterable(data))
        return {key : data}

    @classmethod
    def _remove_comments(cls, content: str) -> str:
        """remove comment lines from the content"""
        content = content.splitlines()
        content = [i for i in content if "COMMENT" not in i ]
        result = "\n".join(content)
        return result

    @classmethod
    def _remove_empty_data(cls, content_dict: str) -> str:
        """remove None s from the parsed content"""
        for k in content_dict:
            content_dict[k] = [i for i in content_dict[k] if i is not None ]

    @classmethod
    def _parse_title(cls, content: str) -> Dict:
        """parse the 'title' section of prmtop file format."""
        result = cls._parse_general(content)
        result["TITLE"] = "".join([i for i in result["TITLE"] if i])
        return result

    @classmethod
    def _parse_atom_name(cls, content: str) -> Dict:
        """parse the 'atom_name' section of prmtop file format."""
        result = cls._parse_general(content)
        cls._remove_empty_data(result)
        return result

    @classmethod
    def _parse_charge(cls, content: str) -> Dict:
        """parse the 'charge' section of prmtop file format.
        NOTE that it is in Amber unit i.e. 18.2223e """
        result = cls._parse_general(content)
        cls._remove_empty_data(result)
        return result

    @classmethod
    def _parse_atomic_number(cls, content: str) -> Dict:
        """parse the 'atomic_number' section of prmtop file format."""
        result = cls._parse_general(content)
        cls._remove_empty_data(result)
        return result

    @classmethod
    def _parse_residue_label(cls, content: str) -> Dict:
        """parse the 'residue_label' section of prmtop file format."""
        result = cls._parse_general(content)
        cls._remove_empty_data(result)
        return result

    @classmethod
    def _parse_residue_pointer(cls, content: str) -> Dict:
        """parse the 'residue_pointer' section of prmtop file format. TODO"""
        result = cls._parse_general(content)
        cls._remove_empty_data(result)
        return result

    @classmethod
    def _parse_residue_number(cls, content: str) -> Dict:
        """parse the 'residue_number' section of prmtop file format."""
        content = cls._remove_comments(content)
        result = cls._parse_general(content)
        cls._remove_empty_data(result)
        return result

    @classmethod
    def _parse_atom_number(cls, content: str) -> Dict:
        """parse the 'atom_number' section of prmtop file format. TODO"""
        content = cls._remove_comments(content)
        result = cls._parse_general(content)
        return result

    @classmethod
    def _parse_atom_element(cls, content: str) -> Dict:
        """parse the 'atom_element' section of prmtop file format. TODO"""
        content = cls._remove_comments(content)
        result = cls._parse_general(content)
        return result

    @classmethod
    def _parse_residue_chainid(cls, content: str) -> Dict:
        """parse the 'residue_chainid' section of prmtop file format. TODO"""
        content = cls._remove_comments(content)
        result = cls._parse_general(content)
        return result

    # region unfullfilled TODO finish these when needed. I didnt make them directly use _parse_general because we dont want to store too much we dont need in the mem.
    @classmethod
    def _parse_atom_bfactor(cls, content: str) -> Dict:
        """parse the 'atom_bfactor' section of prmtop file format. TODO"""
        result = {}

        return result

    @classmethod
    def _parse_atom_occupancy(cls, content: str) -> Dict:
        """parse the 'atom_occupancy' section of prmtop file format. TODO"""
        result = {}

        return result

    @classmethod
    def _parse_pointers(cls, content: str) -> Dict:
        """parse the 'pointers' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_mass(cls, content: str) -> Dict:
        """parse the 'mass' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_atom_type_index(cls, content: str) -> Dict:
        """parse the 'atom_type_index' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_number_excluded_atoms(cls, content: str) -> Dict:
        """parse the 'number_excluded_atoms' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_nonbonded_parm_index(cls, content: str) -> Dict:
        """parse the 'nonbonded_parm_index' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_bond_force_constant(cls, content: str) -> Dict:
        """parse the 'bond_force_constant' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_bond_equil_value(cls, content: str) -> Dict:
        """parse the 'bond_equil_value' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_angle_force_constant(cls, content: str) -> Dict:
        """parse the 'angle_force_constant' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_angle_equil_value(cls, content: str) -> Dict:
        """parse the 'angle_equil_value' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_dihedral_force_constant(cls, content: str) -> Dict:
        """parse the 'dihedral_force_constant' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_dihedral_periodicity(cls, content: str) -> Dict:
        """parse the 'dihedral_periodicity' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_dihedral_phase(cls, content: str) -> Dict:
        """parse the 'dihedral_phase' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_scee_scale_factor(cls, content: str) -> Dict:
        """parse the 'scee_scale_factor' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_scnb_scale_factor(cls, content: str) -> Dict:
        """parse the 'scnb_scale_factor' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_solty(cls, content: str) -> Dict:
        """parse the 'solty' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_lennard_jones_acoef(cls, content: str) -> Dict:
        """parse the 'lennard_jones_acoef' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_lennard_jones_bcoef(cls, content: str) -> Dict:
        """parse the 'lennard_jones_bcoef' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_bonds_inc_hydrogen(cls, content: str) -> Dict:
        """parse the 'bonds_inc_hydrogen' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_bonds_without_hydrogen(cls, content: str) -> Dict:
        """parse the 'bonds_without_hydrogen' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_angles_inc_hydrogen(cls, content: str) -> Dict:
        """parse the 'angles_inc_hydrogen' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_angles_without_hydrogen(cls, content: str) -> Dict:
        """parse the 'angles_without_hydrogen' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_dihedrals_inc_hydrogen(cls, content: str) -> Dict:
        """parse the 'dihedrals_inc_hydrogen' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_dihedrals_without_hydrogen(cls, content: str) -> Dict:
        """parse the 'dihedrals_without_hydrogen' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_excluded_atoms_list(cls, content: str) -> Dict:
        """parse the 'excluded_atoms_list' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_hbond_acoef(cls, content: str) -> Dict:
        """parse the 'hbond_acoef' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_hbond_bcoef(cls, content: str) -> Dict:
        """parse the 'hbond_bcoef' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_hbcut(cls, content: str) -> Dict:
        """parse the 'hbcut' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_amber_atom_type(cls, content: str) -> Dict:
        """parse the 'amber_atom_type' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_tree_chain_classification(cls, content: str) -> Dict:
        """parse the 'tree_chain_classification' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_join_array(cls, content: str) -> Dict:
        """parse the 'join_array' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_irotat(cls, content: str) -> Dict:
        """parse the 'irotat' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_radius_set(cls, content: str) -> Dict:
        """parse the 'radius_set' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_radii(cls, content: str) -> Dict:
        """parse the 'radii' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_screen(cls, content: str) -> Dict:
        """parse the 'screen' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_solvent_pointers(cls, content: str) -> Dict:
        """parse the 'solvent_pointers' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_atoms_per_molecule(cls, content: str) -> Dict:
        """parse the 'atoms_per_molecule' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_box_dimensions(cls, content: str) -> Dict:
        """parse the 'box_dimensions' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_cap_info(cls, content: str) -> Dict:
        """parse the 'cap_info' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_cap_info2(cls, content: str) -> Dict:
        """parse the 'cap_info2' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_polarizability(cls, content: str) -> Dict:
        """parse the 'polarizability' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        result = {}

        return result

    @classmethod
    def _parse_ipol(cls, content: str) -> Dict:
        """parse the 'ipol' section of prmtop file format. TODO"""
        # return cls._parse_general(content)
        # NOTE this is an example that general dont work
        result = {}

        return result
    # endregion
