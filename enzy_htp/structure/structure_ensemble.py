"""This class define one of the core data structure of EnzyHTP: StructureEnsemble.
StructureEnsemble stands for a collection of different geometries of the same enzyme structure or
most commonly, trajectories sampled from a simultion.
As the core data class, StructureEnsemble will be solely designed for **storing, accessing and editing** data.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-28
"""
from typing import List
from .structure import Structure
from .structure_io import StructureParserInterface

class StructureEnsemble:
    """A collection of different geometries of the same enzyme structure.
    Since it is not a rare case that StructureEnsemble can go to >10GB, we will store
    them mostly in files (instead of in memory). And parse them when needed.

    Attributes:
        topology
        top_parser
        coordinate_list
        coord_parser

    Derived properties:
        structures: List[Structure]

    """
    def __init__(self,
                 topology, top_parser: StructureParserInterface,
                 coordinate_list: List, coord_parser: StructureParserInterface) -> None:
        self.topology = topology
        self.top_parser = top_parser
        self.coordinate_list = coordinate_list
        self.coord_parser = coord_parser


    @property
    def structures(self) -> List[Structure]:
        pass
