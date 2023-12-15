"""This class defines the format of MD results in general in EnzyHTP. The goal is to make
MD results engine independent. Efforts may include containing a parsing object from the engine
It is placed here because it should not depend on any knowledge of interface or Structure.

Author: QZ Shao, <shaoqz@icloud.com>
Date: 2023-09-19
"""
from abc import ABC, abstractmethod
from typing import Callable, Dict

from enzy_htp.structure import Structure, StructureEnsemble

class MolDynResult:
    """The unified format of MD result in EnzyHTP.
    The MD will provide 1. a trajectory and 2. some by-product metrics along the
    simulation such as temperature, pressure, etc.
    2 can be calculated from 1 afterwards in most case but are calculated along the way too.
    
    Add more when needed for 1. something that cannot be calculated from traj (could be velocity) or
    2. to avoid redundant heavy calculation.
    
    Attribute:
        traj_file
        traj_parser
        traj_log_file
        traj_log_parser
        last_frame_file # used for next step
        last_frame_parser
    (Note: it is store in file and parser because in most cases you dont want to parse
    a traj because it is too large. Most cases traj uses a lazy parsing scheme that only
    parsed when have to.)
    """
    def __init__(self,
                 traj_file: str,
                 traj_parser: Callable[[str], StructureEnsemble],
                 traj_log_file: str,
                 traj_log_parser: Callable[[str], Dict],
                 last_frame_file: str,
                 last_frame_parser: Callable[[str], Structure],):
        self._traj_file = traj_file
        self._traj_parser = traj_parser
        self._traj_log_file = traj_log_file
        self._traj_log_parser = traj_log_parser
        self._last_frame_file = last_frame_file
        self._last_frame_parser = last_frame_parser

    @property
    def traj_parser(self) -> str:
        return self._traj_file

    @property
    def traj_file(self) -> Callable:
        return self._traj_parser

    @property
    def traj_log_file(self) -> str:
        return self._traj_log_file

    @property
    def traj_log_parser(self) -> Callable:
        return self._traj_log_parser

    @property
    def last_frame_file(self) -> str:
        return self._last_frame_file

    @property
    def last_frame_parser(self) -> Callable:
        return self._last_frame_parser


class MolDynResultEgg(ABC):
    """This class defines the format of md result eggs.
    These eggs are md result place holders before the calculation.
    An example of eggs is file paths, so in Amber's case, it could be
    the .nc file path etc. In general it should be some info to deduce
    where the result will be.
    The motivation of having this class is when running MD through ARMer
    the jobs of non-1st-md-step needs to be configured before previous
    job finishes so that they can know where to find the result of the
    previous step when themselves starts.
    In most cases, the concrete class is a dataclass"""
    pass
