"""Submodule contains code for calculations the RMSD value of a StructureEnsemble instance.
The calculation can be performed with Amber interface.
+ rmsd_calculation()
    Calculate the RMSD value of a StructureEnsemble instance with specified mask pattern.
+ rmsd_of_region()
    Get RMSD value from MD simulation result with a mask region (a list of residues).
+ rmsd_of_structure()
    Get RMSD value from MD simulation result of the whole structure.

Author: Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>

Date: 2024-11-15
"""
# Here put the import lib.
from typing import List, Dict, Union

# Here put EnzyHTP imports.
from enzy_htp import interface as eh_interface
from enzy_htp.preparation import remove_solvent
from enzy_htp.structure import StruSelection
from enzy_htp.structure.structure_ensemble import StructureEnsemble
from enzy_htp.structure.structure_selection import select_stru

def rmsd_calculation(stru_esm: StructureEnsemble, stru_selection: StruSelection) -> List[float]:
    """Calculate the RMSD value of a StructureEnsemble instance with selection.
    
    Args:
        stru_esm (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        stru_selection (StruSelection): A StruSelection instance representing the region for calculating RMSD value.
    """
    return eh_interface.amber.get_rmsd(stru_esm=stru_esm, stru_selection=stru_selection)

def rmsd_of_pattern(stru_esm: StructureEnsemble, region_pattern: str = "all and (not elem H)") -> List[float]:
    """Calculate the RMSD value of a StructureEnsemble instance with specified mask pattern.
    
    Args:
        stru_esm (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        region_pattern (str): A pymol-formatted selection string which defines the region for calculating RMSD value.
    """
    stru = remove_solvent(stru_esm.structure_0, in_place=False)
    stru_sele = select_stru(stru, pattern=region_pattern)
    return rmsd_calculation(stru_esm=stru_esm, stru_selection=stru_sele)

def rmsd_of_structure(stru_esm: StructureEnsemble, include_ligand: bool = True, ca_only: bool = True) -> List[float]:
    """Get RMSD value from MD simulation result of the whole structure. Solvents are not included.
    
    Args:
        stru_esm (StructureEnsemble): A collection of different geometries of the same enzyme structure.
        include_ligand (bool, optional): Indicate if ligands are included in RMSD calculation. Default True.
        ca_only (bool, optional): Indicate if only C-alpha are included in RMSD calculation; otherwise all atoms except hydrogens are included. Default True.
    """
    stru = remove_solvent(stru_esm.structure_0, in_place=False)
    region_pattern = "(not solvent) and (not inorganic)"
    if (not include_ligand):
        ligand_list = list(filter(lambda residue: (residue.is_noncanonical()), stru.residues))
        region_pattern += f"and (not resi {'+'.join([str(ligand.idx) for ligand in ligand_list])})"
    if (ca_only):
        region_pattern += f" and (name CA)"
    else:
        region_pattern += f" and (not elem H)"
    stru_sele = select_stru(stru, pattern=region_pattern)
    return rmsd_calculation(
        stru_esm=stru_esm, 
        stru_selection=stru_sele,
    )
