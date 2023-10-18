"""Here store code that manage NCAA libraries that the default one shared among EnzyHTP.
It is placed under _interface since it is mostly used by interfaces to external software.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-17"""
import re
from typing import List, Tuple, Union

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.structure import (
    Residue,
    PrepinParser,
)

MOL_DESC_FILE_EXT = [".prepin", ".prepi", ".mol2"]
MM_PARM_FILE_EXT = [".frcmod"]

def search_ncaa_parm_file(target_res: Residue, ncaa_lib_path: str) -> Tuple[str, List[str]]: 
    """search for ncaa parm files for {target_res_name} from {ncaa_lib_path}.
    A cache file stores known {file : res_name} relationship in pickle format under the path.
    Args:
        target_res: the target Residue child class instance. (e.g.: Ligand, ModifiedResidue)
    Returns:
        (prepi/mol2_path, [frcmod_path, ...]) if found
        (None, []) if not found"""
    mol_desc_path = None
    frcmod_path_list = []

    # I. initiate the ncaa_lib_mapper {file : res_name} (1-time effort if no new files added)
    ncaa_lib_mapper = {}

    parm_file_list = fs.get_all_file_in_dir(ncaa_lib_path)
    for parm_file in parm_file_list:
        res_name = self.cache_ncaa_lib_mapper.get(parm_file, None) # TODO change this to a cache pickle file instead
        if res_name: # exist in cache
            
            ncaa_lib_mapper[parm_file] = res_name

        else: # not exist in cache - file the res_name
            # == all file support goes here ==
            # prepin
            if fs.get_file_ext(parm_file) in [".prepin", ".prepi"]:
                # 1. find 3-letter name in file
                parm_stru = PrepinParser().get_structure(parm_file)
                res_name = parm_stru.residues[0].name
                if res_name in [None, "UNK"]:
                    # 2. find 3-letter name in filename
                    res_name = self._get_ncaa_parm_file_res_name_from_filename(parm_file)
            # mol2
            elif fs.get_file_ext(parm_file) in [".mol2"]:
                pass
                # # 1. find 3-letter name in file
                # parm_stru = Mol2Parser().get_structure(parm_file)
                # res_name = parm_stru.residues[0].name
                # if res_name in [None, "UNK"]:
                #     # 2. find 3-letter name in filename
                #     res_name = self._get_ncaa_parm_file_res_name_from_filename(parm_file)

            # frcmod*
            elif ".frcmod" in fs.get_file_ext(parm_file):
                # seems frcmod file wont contain residue name in file
                res_name = _get_ncaa_parm_file_res_name_from_filename(parm_file)

            else:
                _LOGGER.warning(
                    f"The file: {parm_file} in ncaa_parm_lib have an unknown extension."
                    "This file will not be considered during parameterization"
                    "currently supported types are: prepin, mol2, frcmod*")
                continue

            # known file type
            if res_name:
                ncaa_lib_mapper[parm_file] = res_name
            else:
                _LOGGER.warning(
                    f"The file: {parm_file} in ncaa_parm_lib does not have an "
                    "associated residue name. This file will not be considered "
                    "during parameterization. Make sure this is what you want! "
                    "You can add residue name by adding it to the corresponding "
                    "part in the file format or use the 3-letter name as filename (in upper case)") 

    self.cache_ncaa_lib_mapper.update(ncaa_lib_mapper) # cache known ones

    # II. assign to target_res
    for file_path, res_name in ncaa_lib_mapper.items():
        if res_name == target_res.name:
            # prepin/mol2/...
            if fs.get_file_ext(file_path) in MOL_DESC_FILE_EXT:
                mol_desc_path = file_path
            # frcmod*/...
            if MM_PARM_FILE_EXT in fs.get_file_ext(file_path):
                frcmod_path_list.append(file_path)

    return mol_desc_path, frcmod_path_list

def _get_ncaa_parm_file_res_name_from_filename(filename: str) -> Union[str, None]:
    """as name desc"""
    base_file_name = fs.base_file_name(filename)
    if re.match("[A-Z][A-Z][A-Z]", base_file_name):
        return base_file_name
    else:
        return None
