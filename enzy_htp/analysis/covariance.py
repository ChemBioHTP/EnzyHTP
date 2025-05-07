"""Submodule contains code for calculations the covariance matrix from a structure ensemble
+ coord_covariance()
    Calculate the atomic coordinate covariance matrix of a StructureEnsemble instance.

Author: Qianzhen Shao <shaoqz@icloud.com>

Date: 2025-05-05
"""
from typing import Dict, Union

from enzy_htp.preparation import remove_solvent
from enzy_htp import interface as eh_interface
from enzy_htp.structure import StruSelection
from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.structure.structure_selection import select_stru
from enzy_htp.core.logger import _LOGGER

def coord_covariance(
        stru_esm: StructureEnsemble, 
        region_pattern: Union[str, tuple] = "polymer and (not elem H)",
        reference_type: str = "average", 
        mass_weighted: bool = False,
        ignore_solvent: bool = True,
        ) -> Dict[str, float]:
    """Calculate the atomic coordinate covariance matrix the region_pattern of a StructureEnsemble
    instance. Covariance is calculated based on the average structure or the first structure.
    Powered by the matrix covar/mwcovar from Cpptraj from AmberTools for now.

    Args:
        stru_esm: 
            A conformational ensemble of a structure.
        region_pattern: 
            A pymol-style selection pattern that defines atoms of the covariance calculation. (i.e., a 3N*3N matrix)
            You can also specify a tuple of two patterns which covariance between atoms described
            in them will be calculated. (i.e., a 3N*3M matrix)
        reference_type:
            control the reference coordinate of calculating the variance of each atom.
            "average" - use an average structure of the ensemble
            "first" - use the first structure of the ensemble
        mass_weigthed:
            control whether the atomic mass is used to weight the covariance.
        ignore_solvent:
            control of solvent are removed before the calculation. True by default to speed up the calculation
            if solvent is the interest of study, set this to False.

    Returns:
        a 2D numpy matrix. i.e., the cartesian coordinate covariance matrix.
            When one pattern is provided it is 3N*3N
            When two pattern is provided it is 3N*3M (N, M are number of atoms from each selection)
    """
    if ignore_solvent:
        stru = remove_solvent(stru_esm.structure_0)
    else:
        stru = stru_esm.structure_0

    if isinstance(region_pattern, tuple): # NOTE the support is already there in AmberInterface just not tested
        _LOGGER.error("specifying two patterns are not supported yet. contact developer to request this feature")
        raise TypeError
    elif isinstance(region_pattern, str):
        stru_sele: StruSelection = select_stru(stru, pattern=region_pattern)
    else:
        _LOGGER.error("only support str or tuple")
        raise TypeError

    return eh_interface.amber.get_coord_covariance(
        stru_esm=stru_esm, 
        stru_selection=stru_sele, 
        reference_type=reference_type, 
        mass_weighted=mass_weighted
    )
