"""

DOI: https://doi.org/10.1021/acs.jpclett.3c02444

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-04-18
"""
from enzy_htp import interface, _LOGGER

from typing import List, Union, Any

from enzy_htp.structure import (Structure, StructureEnsemble, Ligand, PDBParser, Residue)

from plum import dispatch


@dispatch
def spi_metric(
    stru : StructureEnsemble,
    ligand : Ligand,
    pocket_sele: str 
) -> List[float]:
    """
    Args:
        stru:
        ligand:
        pocket_sele:

    Returns:
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
    """
    """

    
    if isinstance(stru, Structure):
        #TODO(CJ): put some kind of thing here
        stru = stru.from_single_stru( stru )

    pocket_sele = " or ".join(map(
        lambda rr: f"( chain {rr.parent.name} and resi {rr.idx} )",
        pocket
    ))
    print(pocket_sele)                

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
    """
    """

    (lig_chain, lig_idx) = ligand.key()
    session = interface.pymol.new_session()
    interface.pymol.load_enzy_htp_stru(session, stru )
    results:List[Any] = interface.pymol.general_cmd(session,[
        ('set', 'dot_solvent', 1),
        ('create', 'ligand', f'chain {lig_chain} and resi {lig_idx}'),
        ('create', 'protein', f'not (chain {lig_chain} and resi {lig_idx})'),
        ('get_area', 'ligand'),
        ('get_area', f'protein and ({pocket_sele})')
    ])
    return results[-2] / results[-1]
