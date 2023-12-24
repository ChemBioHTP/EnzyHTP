"""Here store code that manage NCAA libraries that the default one shared among EnzyHTP.
It is placed under _interface since it is mostly used by interfaces to external software.

The file format in those NCAA libraries should be "RES[_METHOD].ext".
    *RES* is the 3-letter residue name of the NCAA
        for file types like "prepin" that have residue name in the file,
        the value here in filename wont affect anything but it cannot be empty.
    *METHOD* is the method name that is special to the file generation (e.g.: AM1BCC)
        [] means this part is optional
    *.ext* is the extension that is used to judge the file type.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-17"""
import re
import pickle
import time
from pathlib import Path
from typing import List, Tuple, Union

from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp.structure import (
    NonCanonicalBase,
    PrepinParser,
    Mol2Parser,
)

MOL_DESC_FILE_EXT = [".prepin", ".prepi", ".mol2"]
"""extensions of molecule description file types"""

MM_PARM_FILE_EXT = [".frcmod"]
"""extensions of parameter file types"""

PARM_METHOD_LIST = ["AM1BCC-GAFF2", "RESP-GAFF2", "AM1BCC-GAFF", "RESP-GAFF", "MCPB", "any"]
"""a list of keywords for supported parm_method as a search target. If 'any' is used, any method
will be matched including None or ''. * GAFF is used if not specificed."""

def search_ncaa_parm_file(target_res: NonCanonicalBase, target_method: str,
                          ncaa_lib_path: str, use_cache: bool= True) -> Tuple[str, List[str]]:
    """search for ncaa parm files for {target_res_name} with {target_method} from {ncaa_lib_path}.
    A cache file will be made to store known {file : (res_name, method)} relationship in pickle
    format under the path.
    Args:
        target_res:
            the target Residue child class instance. (e.g.: Ligand, ModifiedResidue)
        target_method:
            the desired method for the parm file
            (if "any" is supplied, any method is accepted, first found is used)
    Returns:
        (prepi/mol2_path, [frcmod_path, ...]) if found
        (None, []) if not found"""
    mol_desc_path = None
    frcmod_path_list = []
    # 0. san check
    if target_method not in PARM_METHOD_LIST:
        _LOGGER.error(f"target 'method' ({target_method}) of search is not valid."
                      f"Supported list: {PARM_METHOD_LIST}."
                      "If you are not searching for a specific method, please use 'any'."
                      "To register it as a supported method - add it to ncaa_library.py::PARM_METHOD_LIST.")
        raise ValueError
    if not Path(ncaa_lib_path).exists():
        _LOGGER.error(f"supplied ncaa_lib_path does not exist: {ncaa_lib_path}")
        raise ValueError

    # I. initiate a ncaa_lib_mapper {file : (res_name, method)} (1-time effort if no new files added)
    ncaa_lib_mapper = {}
    # 1. init cache from file
    cache_file_path = f"{ncaa_lib_path}/.cache_lib_mapper.pickle"
    if not Path(cache_file_path).exists():
        with open(cache_file_path, "wb") as of:
            pickle.dump({}, of)
    with open(cache_file_path, "rb") as f:
        while fs.is_locked(f): # wait if the file is writing by other workflow
            time.sleep(0.1)
        cache_ncaa_lib_mapper = pickle.load(f)

    # 2. deduce (res_name, method) for each file
    parm_file_list = fs.all_file_in_dir(ncaa_lib_path)
    for parm_file in parm_file_list:

        res_name, parm_method = cache_ncaa_lib_mapper.get(parm_file, (None, None))

        if res_name and use_cache: # exist in cache (because we allow parm_method to be None)
            _LOGGER.debug(f"using cache for {parm_file} : {(res_name, parm_method)} ")
            ncaa_lib_mapper[parm_file] = (res_name, parm_method)

        else: # not exist in cache - deduce (res_name, method)

            # == all file support goes here ==
            # prepin
            if fs.get_file_ext(parm_file) in [".prepin", ".prepi"]:
                # 2.1 find 3-letter name in file
                parm_stru = PrepinParser().get_structure(parm_file)
                res_name = parm_stru.residues[0].name
                if res_name in [None, "UNK"]:
                    # 2.2 find 3-letter name in filename
                    res_name = _get_res_name_from_filename(parm_file)
                # 2.3 find method name in filename
                parm_method = _get_parm_method_from_filename(parm_file)

            # mol2
            elif fs.get_file_ext(parm_file) in [".mol2"]:
                # 2.1. find 3-letter name in file
                parm_stru = Mol2Parser().get_structure(parm_file)
                res_name = parm_stru.residues[0].name
                if res_name in [None, "UNK"]:
                    # 2.2 find 3-letter name in filename
                    res_name = _get_res_name_from_filename(parm_file)
                # 2.3 find method name in filename
                parm_method = _get_parm_method_from_filename(parm_file)

            # frcmod*
            elif ".frcmod" in fs.get_file_ext(parm_file):
                # seems frcmod file wont contain residue name in file
                res_name = _get_res_name_from_filename(parm_file)
                parm_method = _get_parm_method_from_filename(parm_file)
            else:
                _LOGGER.warning(
                    f"The file: {parm_file} in ncaa_parm_lib have an unknown extension."
                    "This file will not be considered during parameterization"
                    "currently supported types are: prepin, mol2, frcmod*")
                continue

            # known file type
            if res_name:
                ncaa_lib_mapper[parm_file] = (res_name, parm_method)
            else:
                _LOGGER.warning(
                    f"The file: {parm_file} in ncaa_parm_lib does not have an "
                    "associated residue name. This file will not be considered "
                    "during parameterization. Make sure this is what you want! "
                    "You can add residue name by adding it to the corresponding "
                    "part in the file format or use the 3-letter name as filename (in upper case)")
            if not parm_method:
                _LOGGER.debug(
                    f"The file: {parm_file} in ncaa_parm_lib does not have an "
                    "associated generation method. This file *will be* considered but"
                    "will only be matched when target_method is 'any'"
                    " Make sure this is what you want!"
                    "You can add 'method' by adding it to the corresponding name"
                    " in the filename (see ncaa_library.py::PARM_METHOD_LIST)")

    # 3. update cache (only when there are sth. new/different)
    if cache_ncaa_lib_mapper != ncaa_lib_mapper:
        cache_ncaa_lib_mapper.update(ncaa_lib_mapper)
        with open(cache_file_path, "wb") as of:
            fs.lock(of)
            pickle.dump(cache_ncaa_lib_mapper, of)
            fs.unlock(of)

    # II. assign to (target_res, target_method)
    for file_path, (res_name, parm_method) in ncaa_lib_mapper.items():
        if res_name == target_res.name and (target_method in [parm_method,"any"]):
            # prepin/mol2/...
            if fs.get_file_ext(file_path) in MOL_DESC_FILE_EXT:
                mol_desc_path = file_path
            # frcmod*/...
            for parm_ext in MM_PARM_FILE_EXT:
                if parm_ext in fs.get_file_ext(file_path):
                    frcmod_path_list.append(file_path)
                    break

    return mol_desc_path, frcmod_path_list

def _get_res_name_from_filename(filename: str) -> Union[str, None]:
    """as name desc"""
    base_file_name = fs.base_file_name(filename)
    res_name_part = base_file_name.split("_")[0]
    if re.match("[A-Z0-9][A-Z0-9][A-Z0-9]", res_name_part):
        return res_name_part
    else:
        return None

def _get_parm_method_from_filename(filename: str) -> Union[str, None]:
    """as name desc"""
    base_file_name_parts = fs.base_file_name(filename).split("_")
    if len(base_file_name_parts) < 2:
        return None
    method = base_file_name_parts[1]
    if method not in PARM_METHOD_LIST:
        _LOGGER.warning(f"NCAA lib file 'method' ({method}) does not match any supported ones. ({PARM_METHOD_LIST}) "
                        "It will only match 'any' in the search. Make sure this is what you want!")
    return method
