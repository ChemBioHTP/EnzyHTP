"""Generation/construction of Structure objects from XYZ files and exporting it vice versa

Author: Sebastian Stull, <sebastian.l.stull@vanderbilt.edu>
Date: 2024-12-14
"""
import os
from typing import List, Union

from plum import dispatch

import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_region.structure_region import StructureRegion
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
    def get_file_str(self, stru: Union[Structure, StructureRegion]) -> str:
        result = f"{str(len(stru.atoms))}{os.linesep}{os.linesep}"
        for at in stru.atoms:
            result += f"{at.element} {at.coord[0]} {at.coord[1]} {at.coord[2]}{os.linesep}"
        result  = result[:-1]
        return result