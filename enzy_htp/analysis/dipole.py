"""Submodule contains code for calculating the molecular dipole in an enzyme system.
Current supported APIs:
+ bond_dipole()
    calculate the bond dipole based on the electronic_structure/wavefunction

Author: QZ Shao <shaoqz@icloud.com>

Date: 2024-01-12
"""
from typing import Tuple, Dict
import numpy as np

from enzy_htp.core.logger import _LOGGER
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.structure import Atom
from enzy_htp import interface

def bond_dipole(
        ele_stru: ElectronicStructure,
        atom_1: Atom,
        atom_2: Atom,
        method: str = "LMO-Multiwfn",
        # execution config
        work_dir: str="./bond_dipole",
        keep_in_file: bool=False,
        # armer config
        cluster_job_config: Dict = None,
        job_check_period: int= 30, # s
        **kwargs,
        ) -> Tuple[float, np.array]:
    """calculate bond dipole based on an ElectronicStructure.
    (TODO in the case of using AMOEBA, probably also stru instead?)
    
    Args:
        ele_stru:
            the ElectronicStructure of a QM region containing the
            target bond.
        atom_1:
            the 1st atom defining target bond (in the region)
        atom_2:
            the 2st atom defining target bond (in the region)
        method:
            the method keyword specifying the algorithm & software
            of the dipole calculation. (see Details for supported)
        work_dir:
            the working dir that contains all the files in the SPE process
        keep_in_file:
            whether keep the input file of the calculation
        cluster_job_config: 
            the config for cluster_job. If None is used, the calculation
            will be run locally.
        job_check_period:
            the time cycle for update job state change (Unit: s)

        (see more options in detailed engines)
    
    Returns:
        (dipole_norm_signed, dipole_vec)
            *dipole_norm_signed*
                is the signed norm of the dipole according to its projection
                to be bond vector.
            *dipole_vec*
                is the vector of the dipole

            Dipole positive Direction: negative(-) to positive(+).
            Result direction: atom_1 -> atom_2

    Details:
        Supported methods:
            "LMO-Multiwfn"
            > Multiwfn manual 3.22/4.19.4
            > Lu, T.; Chen, F., Multiwfn: A multifunctional wavefunction analyzer. J. Comput. Chem. 2012, 33 (5), 580-592.
            2-center LMO dipole is defined by the deviation of the eletronic mass center
            relative to the bond center."""
    # san check
    if not isinstance(ele_stru, ElectronicStructure):
        _LOGGER.error(f"ele_stru can only be an ElectronicStructure(). found: {ele_stru}")
        raise TypeError
    if method not in BOND_DIPOLE_METHODS:
        _LOGGER.error(f"method ({method}) not supported. Supported: {BOND_DIPOLE_METHODS.keys()}")
        raise ValueError

    result = BOND_DIPOLE_METHODS[method](
        ele_stru, atom_1, atom_2, work_dir, keep_in_file,
        cluster_job_config = cluster_job_config, 
        job_check_period = job_check_period,
        **kwargs)

    return result

BOND_DIPOLE_METHODS = {"LMO-Multiwfn": interface.multiwfn.get_bond_dipole}
"""the method keyword mapping for bond dipole calculations."""
