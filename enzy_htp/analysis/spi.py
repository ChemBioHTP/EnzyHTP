"""Calculates the substrate positioning index (spi or spi metric) which uses solvent accessible surface area (SASA)
to calculate a ratio between the ligand and defined pocket. In other words, spi = SASA_ligand / SASA_pocket. Functionality 
is accessed through spi_metric(), which is  overloaded with the below signatures:

    + spi_metric( StructureEnsemble, str, str ) -> List[float]
    + spi_metric( StructureEnsemble, Ligand, List[Residue] ) -> List[float]
    + spi_metric( StructureEnsemble, Ligand, str) -> List[float]
    + spi_metric( Structure, Ligand, str) -> float

The most used API will typically be spi_metric( StructureEnsemble, Ligand, List[Residue] ). To get the List[Residue], consider
using enzy_htp.structure.structure_operation.closest_n_residues().

If the ligand consists of multiple Residue instances, consider using spi_metric( StructureEnsemble, str, str ), 
which enables you to define your ligand with pymol-formatted selection pattern.

Note that all SASA calculations are done using pymol's python API, using the function PyMolInterface.get_spi().
Citation for spi: DOI:https://doi.org/10.1021/acs.jpclett.3c02444. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>; Zhong, Yinjie <yinjie.zhong@vanderbilt.edu>
Date: 2025-01-28
"""
from enzy_htp import interface, _LOGGER

from typing import List, Union, Any

from enzy_htp.structure import (Structure, StructureEnsemble, Ligand, Residue)

from plum import dispatch


@dispatch
def spi_metric(
    stru: StructureEnsemble,
    ligand_pattern: str,
    pocket_pattern: str 
) -> List[float]:
    """Calculates the spi metric for a StructureEnsemble using a pymol-formatted pocket selection pattern str. 
    Note: The pocket_pattern expression is applied to each Structure instances in the StructureEnsemble instance,
    so it is recommended to use a selection which specifies individual residues, 
    as the expression defined based on distance (to the ligand) will be re-evaluated for each Structure instance's geometry.

    Args:
        stru (StructureEnsemble): The `StructureEnsemble` instance to analyze.
        ligand_pattern (str): A pymol-formatted selection pattern str which defines the ligand (the numerator) in the spi metric.
        pocket_pattern (str): A pymol-formatted selection pattern str which defines the pocket (the denominator) in the spi metric.

    Returns:
        A List[float] have the same length of `StructureEnsemble.structures`, with each value being a different spi metric for each `Structure` instance.
    """
    result = list()
    for single_stru in stru.structures():
        result.append(
            interface.pymol.get_spi(single_stru, ligand_pattern, pocket_pattern)
        )
    return result

@dispatch
def spi_metric(
    stru : StructureEnsemble,
    ligand : Ligand,
    pocket_pattern: str 
) -> List[float]:
    """Calculates the spi metric for a StructureEnsemble using a pymol-formatted pocket sele str. 
    Note: The pocket_pattern expression is applied to each Structure(), so it is recommended to use a selection which specifies
    individual Residue()'s as the expression defined based on distance will be re-evaluated for each Structure() geometry.

    Args:
        stru (StructureEnsemble): The StructureEnsemble() to analyze.
        ligand (Ligand): The Ligand() which will be the numerator in the spi metric.
        pocket_pattern (str): A pymol-formatted sele str which defines the denominator in the spi metric.

    Returns:
        A List[float] of len(StructureEnsemble.structures), with each value being a different spi metric for the nth Structure(). 
    """
    
    if isinstance(stru, Structure):
        #TODO(CJ): put some kind of thing here
        stru = stru.from_single_stru( stru )

    result = list()
    for ss in stru.structures():
        result.append(
            spi_metric( ss, ligand, pocket_pattern )
        )

    return result

@dispatch
def spi_metric(
    stru : StructureEnsemble,
    ligand : Ligand,
    pocket: List[Residue],
) -> List[float]:
    """Calculates the spi metric for a StructureEnsemble using a List[Residue] that define the denominator of the 
    metric. This is the PREFERRED method for calculating spi.

    Args:
        stru (StructureEnsemble): The StructureEnsemble() to analyze.
        ligand (Ligand): The Ligand() which will be the numerator in the spi metric.
        pocket (List[Residue]): A List[Residue] which serves as the pocket or SASA denominator of the spi metric.

    Returns:
        A List[float] of len(StructureEnsemble.structures), with each value being a different spi metric for the nth Structure(). 
    """
    pocket_sele = " or ".join(map(
        lambda rr: f"( chain {rr.parent.name} and resi {rr.idx} )",
        pocket
    ))

    result = list()
    for ss in stru.structures():
        result.append(
            spi_metric(ss, ligand, pocket_sele)
        )

    return result

@dispatch
def spi_metric( stru: Structure, 
            ligand: Ligand,
            pocket_sele:str ) -> float:
    """Calculates the spi metric for a single Structure using a pymol-formatted pocket sele str. 
    Note: The pocket_sele expression is applied to the Structure(), so it is recommended to use a selection which specifies
    individual Residue()'s. Relies on the PyMolInterfaace.get_spi() method.

    Args:
        stru: The Structure() to analyze.
        ligand: The Ligand() which will be the numerator in the spi metric.
        pocket_sele: A pymol-formatted sele str which defines the denominator in the spi metric.

    Returns:
        A List[float] of len(StructureEnsemble.structures), with each value being a different spi metric for the nth Structure(). 
    """
   
    return interface.pymol.get_spi(stru, ligand, pocket_sele)


