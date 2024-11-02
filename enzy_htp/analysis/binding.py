"""Submodule contains code for calculations the binding energy of a ligand in the enzyme.
The energy can be either enthalpy or Gibbs free energy depending on the algorithm used
+ binding_energy()
    calculate the binding energy of a ligand in a protein from a Structure or an ensemble
    of Structure()s

Author: QZ Shao <shaoqz@icloud.com>

Date: 2024-05-30
"""
from typing import List, Dict, Tuple, Union, Callable
from functools import partial

from enzy_htp import interface
from enzy_htp.structure.structure_selection import select_stru
from enzy_htp.structure import Structure, Ligand, StructureEnsemble
from enzy_htp.core import _LOGGER
from enzy_htp.core.job_manager import ClusterJob

def binding_energy(
        # Input Data
        stru: Union[Structure, StructureEnsemble],
        ligand: str,
        # Config
        method: str = "mmpbsa_amber",
        cluster_job_config: Dict= None,
        job_check_period: int= 30, # s
        work_dir: str="./binding",
        keep_in_file: bool=False,
        **kwargs,
        ) -> List[float]:
    """Calculate the binding energy of {ligand} in {stru}

    Args:
        stru:
            the target protein as the receptor of the calculation represented as Structure()
        ligand:
            the target ligand of the calculation represented as a selection pattern.
            Note that the ligand has to be part of Structure().
            Note that the ligand can be a small molecule or a protein   
        method:
            the algorithm of the binding energy calculation. (see Details)
        cluster_job_config: 
            the config for cluster_job. (also enables running via a ClusterJob/ARMer)
        job_check_period:
            the time cycle for update job state change (Unit: s)
        work_dir:
            the working dir that contains all the files in the process
        keep_in_file:
            whether keep the input file of the calculation
    (when method="mmpbsa_amber")
        igb:
        TODO

    Returns:
        a list of float for the value of the binding energy
        (if stru is Structure() only 1 value will be in the list)

    Details:
        The binding energy assesses how favorable would a ligand bind into a protein. It is the energy 
        difference between the bound state (protein-ligand) and the unbound state.

    Available implemtations in this field:
        (supported as a "method" keyword)
        - MMPB/GBSA-1A (https://www.tandfonline.com/doi/full/10.1517/17460441.2015.1032936, https://ambermd.org/tutorials/advanced/tutorial3/index.php, https://pubs.acs.org/doi/10.1021/ar000033j)
            keyword: "mmpbsa_amber" or "mmgbsa_amber"
            calculate the MMPB/GBSA binding energy from a single trajectory/StructureEnsemble of the complex using AmberTools
            For each frame from the trajectory, 
            MM:
                the MM energy of the dry complex, dry ligand, and dry protein is calculated 
                and the dry binding energy is calculated by E_dry = E_solvated_complex - E_dry_complex - E_dry_ligand
            PB/GB:
                implicit solvent model is used to calculated the solvation energy for each component
            SA:
                the surface area approximation is used for calculating the hydrophobic component of the solvation free energy
            The total binding energy in solvent is E_dry + E_solv_complex - E_solv_ligand - E_solv_protein
            * Most MMPB/GBSA calculation result approximates the binding enthalpy
            * MMPB/GBSA-1A is bad at calculeting protein-protein interaction (https://pubs.rsc.org/en/content/articlelanding/2016/cp/c6cp03670h)

        (not supported yet)
        - Rosetta Binding Energy
        - Free Energy Perturbation
            + slow
            + accurate
        - AI based method
            AP-net (support under development)"""
 
    # handle the stru dispatch
    if isinstance(stru, Structure):
        # TODO convert to a single frame StructureEnsemble
        raise Exception("TODO. contact the developer if you need this")
    elif isinstance(stru, StructureEnsemble):
        stru_esm = stru
    else:
        _LOGGER.error(f"stru can only be an Structure() or StructureEnsemble(). found: {stru}")
        raise TypeError

    if method not in BINDING_ENERGY_METHODS:
        _LOGGER.error(f"method ({method}) not supported. Supported: {BINDING_ENERGY_METHODS.keys()}")
        raise ValueError

    ligand = select_stru(stru_esm.structure_0, ligand) # the selection has to be resolved here due to the import-layer structure (otherwise casuing loop importing)

    result = BINDING_ENERGY_METHODS[method](
        stru_esm, ligand, 
        keep_in_file = keep_in_file,
        work_dir = work_dir,
        cluster_job_config = cluster_job_config, 
        job_check_period = job_check_period,
        **kwargs)

    return result

BINDING_ENERGY_METHODS: Dict[str, Callable] = {
    "mmpbsa_amber" : partial(interface.amber.get_mmpbgbsa_energy, solvent_model="pbsa"),
    "mmgbsa_amber" : partial(interface.amber.get_mmpbgbsa_energy, solvent_model="gbsa"),
}
