"""Generation/construction of Structure() objects from .prepin files and exporting these objects to this file format. 
Definition of .prepin format (https://ambermd.org/doc/prep.html). All parsing is done within enzy_htp using 
this parser only. The PrepinParser has no private data and serves as a namespace for .prepin I/O conversion functions.

Author: Qianzhen Shao <shaoqz@icloud.com>
Date: 2023-10-17
"""
from ._interface import StructureParserInterface
from ..structure import Structure

class PrepinParser(StructureParserInterface):
    """the parser for AmberMD prepin files"""

    def __init__(self) -> None:  # pylint: disable=super-init-not-called
        """pass"""
        pass

    @classmethod
    def get_structure(cls, path: str) -> Structure:
        """
        Converting a .prepin file (as its path) into the Structure()
        Arg:
            path:
                the file path of the PDB file
        Return:
            Structure()
        """

