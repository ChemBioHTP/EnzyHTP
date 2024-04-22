"""Calculates the substrate positioning index (spi or spi metric) which uses solvent accessible surface area (SASA)
to calculate a ratio between the ligand and defined pocket. In other words, spi = SASA_ligand / SASA_pocket. Functionality 
is accessed through spi_metric(), which is  overloaded with the below signatures:

    + spi_metric( StructureEnsemble, Ligand, List[Residue] ) -> List[float]
    + spi_metric( StructureEnsemble, Ligand, str) -> List[float]
    + spi_metric( Structure, Ligand, str) -> float

The most used API will typically be spi_metric( StructureEnsemble, Ligand, List[Residue] ). To get the List[Residue], consider
using enzy_htp.structure.structure_operation.closest_n_residues(). Note that all SASA calculations are done using pymol's
python API. Citation for spi: DOI:https://doi.org/10.1021/acs.jpclett.3c02444.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-04-18
"""
from enzy_htp import interface, _LOGGER

from typing import List, Union, Any

from enzy_htp.structure import (Structure, StructureEnsemble, Ligand, PDBParser, Residue)

from plum import dispatch
from enzy_htp._interface.pymol_interface import OpenPyMolSession


@dispatch
def spi_metric(
    stru : StructureEnsemble,
    ligand : Ligand,
    pocket_sele: str 
) -> List[float]:
    """Calculates the spi metric for a StructureEnsemble using a pymol-formatted pocket sele str. 
    Note: The pocket_sele expression is applied to each Structure(), so it is recommended to use a selection which specifies
    individual Residue()'s as the expression will be re-evaluated for each Structure() geometry.

    Args:
        stru: The StructureEnsemble() to analyze.
        ligand: The Ligand() which will be the numerator in the spi metric.
        pocket_sele: A pymol-formatted sele str which defines the denominator in the spi metric.

    Returns:
        A List[float] of len(StructureEnsemble.structures), with each value being a different spi metric for the nth Structure(). 
    """
    
    if isinstance(stru, Structure):
        #TODO(CJ): put some kind of thing here
        stru = stru.from_single_stru( stru )

    result = list()
    for ss in stru.structures:
        result.append(
            spi_metric( ss, ligand, pocket_sele )
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
        stru: The StructureEnsemble() to analyze.
        ligand: The Ligand() which will be the numerator in the spi metric.
        pocket: A List[Residue] which serves as the pocket or SASA denominator of the spi metric.

    Returns:
        A List[float] of len(StructureEnsemble.structures), with each value being a different spi metric for the nth Structure(). 
    """
    pocket_sele = " or ".join(map(
        lambda rr: f"( chain {rr.parent.name} and resi {rr.idx} )",
        pocket
    ))

    result = list()
    for ss in stru.structures:
        result.append(
            spi_metric( ss, ligand, pocket_sele )
        )

    return result

@dispatch
def spi_metric( stru: Structure, 
            ligand: Ligand,
            pocket_sele:str ) -> float:
    """Calculates the spi metric for a single Structure using a pymol-formatted pocket sele str. 
    Note: The pocket_sele expression is applied to the Structure(), so it is recommended to use a selection which specifies
    individual Residue()'s.

    Args:
        stru: The Structure() to analyze.
        ligand: The Ligand() which will be the numerator in the spi metric.
        pocket_sele: A pymol-formatted sele str which defines the denominator in the spi metric.

    Returns:
        A List[float] of len(StructureEnsemble.structures), with each value being a different spi metric for the nth Structure(). 
    """
    
    (lig_chain, lig_idx) = ligand.key()
    pi = interface.pymol
    with OpenPyMolSession(pi) as pms:
        interface.pymol.load_enzy_htp_stru(session, pms )
        results:List[Any] = interface.pymol.general_cmd(pms,[
            ('set', 'dot_solvent', 1),
            ('create', 'ligand', f'chain {lig_chain} and resi {lig_idx}'),
            ('create', 'protein', f'not (chain {lig_chain} and resi {lig_idx})'),
            ('get_area', 'ligand'),
            ('get_area', f'protein and ({pocket_sele})')
        ])
        return results[-2] / results[-1]
