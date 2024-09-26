"""This class define one of the core data structure of EnzyHTP: StructureEnsemble.
StructureEnsemble stands for a collection of different geometries of the same enzyme structure or
most commonly, trajectories sampled from a simultion.
As the core data class, StructureEnsemble will be solely designed for **storing, accessing and editing** data.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-28
"""
from __future__ import annotations
from typing import List, Generator, Callable, Tuple
from copy import deepcopy

from .structure import Structure
from .structure_io import StructureParserInterface
from . import structure_operation as stru_oper
from enzy_htp.core.general import get_itself
# from enzy_htp import interface

from pandas import DataFrame, concat
from numpy import ndarray, asarray

# amber_interface = interface.amber

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

    def __init__(
        self,
        topology: str, 
        top_parser: Callable[[str], Structure], 
        coordinate_list: str,
        coord_parser: Callable[[str], Generator[List[Tuple[float]], None, None]]
    ) -> None:
        self._topology = topology
        self.top_parser = top_parser
        self.coordinate_list = coordinate_list
        self.coord_parser = coord_parser

    def structures(self, remove_solvent: bool=False) -> Generator[Structure]:
        """get a Generator of all geometries in the ensemble
        as Structure()s"""
        stru = deepcopy(self.topology)
        if remove_solvent:
            stru_oper.remove_solvent(stru)
            stru_oper.remove_counterions(stru)

        for this_coord in self.coord_parser(
                self.coordinate_list,
                remove_solvent=remove_solvent
            ):
            result = deepcopy(stru)
            if remove_solvent:
                stru_oper.remove_solvent(result)
                stru_oper.remove_counterions(result)
            result.apply_geom(this_coord)
            yield result

    @property
    def topology(self) -> Structure:
        """getter for topology"""
        if isinstance(self._topology, Structure):
            return self._topology
        else:
            return self.top_parser(self._topology)

    @property
    def structure_0(self) -> Structure:
        """getter for the 1st Structure (state) in the ensemble."""
        coord_0 = next(self.coord_parser(self.coordinate_list))
        result = deepcopy(self.topology)
        result.apply_geom(coord_0)
        return result
    
    # def average(self, auto_image=True) -> Structure:
    #     """A structure that represents the average coordinates of the entire Structure Ensemble.
        
    #     Args:
    #         remove_solvent (bool, optional): Indicate if solvents in the structure ensemble are to be removed.
    #     """
        # result = deepcopy(self.topology)
        # if (remove_solvent):
        #     stru_oper.remove_solvent(result)
        #     stru_oper.remove_counterions(result)

        # atom_key_list = [atom.key for atom in result.atoms]
        # coord_df = DataFrame(columns=atom_key_list)

        # structure_generator = self.structures(remove_solvent=remove_solvent)
        # for structure_state in structure_generator:
        #     atom_coord_dict = dict()
        #     for atom in structure_state.atoms:
        #         atom_coord_dict[atom.key] = asarray(atom.coord)
        #         continue
        #     coord_df = concat([coord_df, DataFrame([atom_coord_dict])], ignore_index=True)
        #     continue
        
        # coord_series = coord_df.sum() / len(coord_df)
        # coord_list = [coord.tolist() for coord in coord_series]
        # result.apply_geom(coord_list)
        # return result
        
        # average_structure = amber_interface.get_average_structure(traj_path=self.coordinate_list, prmtop_path=self._topology, auto_image=auto_image)
        # return average_structure
        # pass
    
    @classmethod
    def from_single_stru(cls, stru: Structure) -> StructureEnsemble:
        """create an ensemble of 1 snapshot from a Structure instance"""
        return cls(
            topology=stru,
            top_parser=get_itself,
            coordinate_list=[stru],
            coord_parser=iter,
        )

    # region == special ==
    def __iter__(self):
        return self.structures()
    # endregion
