"""Generation/construction of Structure objects from PDB files and exporting it vice versa
Definition of PDB file format (http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)
Also contains util function that requires knoweledge of PDB format

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2022-08-01
"""
from typing import List, Union

import enzy_htp.core.file_system as fs
from ._interface import StructureParserInterface
from ..atom import Atom
from ..structure import Structure


class XYZParser(StructureParserInterface):

    def __init__(self) -> None:  # pylint: disable=super-init-not-called
        """pass"""
        pass

    @classmethod
    def get_structure(self):
        pass    
    
    @classmethod
    def get_structure(self, path: str) -> Union[Structure, tuple]:
        pass

    @classmethod
    def get_file_str(self, stru: Structure) -> str:
        pass

    @classmethod
    def make_xyz_from_atoms(self, file_path: str, atoms: List[Atom]):
        lines = [
            str(len(atoms)),
            "",]
        for at in atoms:
            lines.append(f"{at.element} {at.coord[0]} {at.coord[1]} {at.coord[2]}")
        fs.write_lines(file_path, lines)
        return file_path