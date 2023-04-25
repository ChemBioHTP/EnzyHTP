"""Opertions that checks Structure().

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-04-24
"""
from enzy_htp.core.logger import _LOGGER
from ..structure import Structure, Solvent, Chain, Residue, Atom

def check_topology_error(stru: Structure):
    """check {stru} for topology error.
    i.e.: rings in structure should not be circling on other bonds.
    An example of this error is in https://github.com/ChemBioHTP/EnzyHTP/issues/110"""
    print(stru) #TODO(eod)

