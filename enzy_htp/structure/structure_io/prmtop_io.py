"""Generation/construction of Structure() objects from .prmtop files and exporting these objects to this file format. 
Definition of .prmtop format (http://ambermd.org/FileFormats.php#topology). All parsing is done within enzy_htp using 
this parser only. The PrmtopParser has no private data and serves as a namespace for .prmtop I/O conversion functions.
TODO

Author: Qianzhen Shao <shaoqz@icloud.com>
Date: 2023-10-28
"""
from typing import Dict, List

from ._interface import StructureParserInterface
from ..structure import Structure, convert_res_to_structure

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import FileFormatError

class PrmtopParser(StructureParserInterface):
    """the parser for Gaussian prmtop files"""

    def __init__(self) -> None:  # pylint: disable=super-init-not-called
        """pass"""
        pass

    @classmethod
    def get_structure(cls, path: str) -> Structure:
        """Converting a .prmtop file (as its path) into the Structure()
        Arg:
            path:
                the file path of the prmtop file
        Return:
            Structure()"""
        pass # TODO: add when need

    @classmethod
    def get_file_str(cls, stru: Structure) -> str:
        """convert a Structure() to .prmtop file content. Only 1 residue unit is allowed in the stru"""
        pass # TODO: add when need

    @classmethod
    def save_structure(cls, out_path: str, stru: Structure) -> str:
        """save a Structure() to .prmtop file. Only 1 residue unit is allowed in the stru.
        return the out_path"""
        pass

    @classmethod
    def _parse_prmtop_file(cls, path: str) -> Dict:
        """parse prmtop file to a data dictionary."""
        pass # TODO: add when need

    