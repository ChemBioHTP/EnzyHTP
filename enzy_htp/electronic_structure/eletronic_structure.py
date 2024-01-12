"""electronic structure module of EnzyHTP. Defines the EletronicStructure
for electronic structure/wavefunctions of a molecule decribed by Structure().

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-12-28"""
import os
from typing import List, Any

from enzy_htp.structure import Structure

class EletronicStructure:
    """This class defines electronic structure/wavefunctions of a 
    molecule decribed by Structure() in EnzyHTP. It is stored in a
    form of files as the number and size of them could go very large.
    The key properties such as MO coefficients will be deduce in a lazy
    determining + cache manner.
    Attributes:
        energy_0: float
            ground state energy
        geometry: Structure
            geometry
        mo: str
            the file path/externel object of MO
        mo_parser: Any
            the adaptor that covert mo to desired properities

    Lazy Attributes:
        mo_coeff: List[List[float]]
            the coefficient matrix of MOs
        mo_occ: List[float]
            the occupency of MOs
        basis_set: str
            the basis function set
    
    Properities:
        TODO"""
    def __init__(self,
        energy_0: float,
        geometry: Structure,
        mo: str,
        mo_parser: Any,
        source: str = None,):
        self._energy_0 = energy_0
        self._geometry = geometry
        self._mo = mo
        self._mo_parser = mo_parser
        self._source = source

    @property
    def energy_0(self) -> float:
        """ground state energy"""
        return self._energy_0

    @property
    def geometry(self) -> Structure:
        """geometry"""
        return self._geometry

    @property
    def mo(self) -> str:
        """the MO file"""
        return self._mo

    
    @property
    def mo_coeff(self) -> List[List[float]]:
        """the coefficient matrix of MOs"""
        return self._mo_parser.get_mo_coeff(self._mo)
    
    @property
    def mo_occ(self) -> List[float]:
        """the occupency of MOs"""
        return self._mo_parser.get_mo_occ(self._mo)
    
    @property
    def basis_set(self) -> str:
        """the basis function set"""
        return self._mo_parser.get_basis_set(self._mo)
    
    @property
    def source(self) -> str:
        """the source of this result"""
        return self._source
    
    # region == Special ==
    def __str__(self) -> str:
        out_line = [
            f"<ElectronicStructure object at {hex(id(self))}>",
        ]
        out_line.append("ElectronicStructure(")
        out_line.append(f"    energy_0: {self.energy_0}")
        out_line.append(f"    MO source: {self.mo}")
        out_line.append(f"    MO parser: {self._mo_parser}")
        out_line.append(")")
        return os.linesep.join(out_line)

    # endregion
