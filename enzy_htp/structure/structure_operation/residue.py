"""Operations on Residue(). Functions in this module take Residue()'s as input and perform 
operations on them.

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-03-20
"""

from typing import Tuple, Union, List

from enzy_htp.core.logger import _LOGGER
from enzy_htp.core import math_helper as mh
import enzy_htp.chemical as chem
from ..structure import Structure, Residue, Atom
from .connectivity import init_connectivity


def deprotonate_residue(residue: Residue, target_atom: Union[None, Atom] = None) -> None:
    """
    Deprotonates the {residue} on the {target_atom} (if provided).
    Removes the acidic hydrogen attached to the {target_atom} and changes the residue name
    corresponding to chem.residue.DEPROTONATION_MAPPER or /resource/ProtonationState.cdx. Performs
    changes in place.
    Args:
        residue: The Residue() to deprotonate.
        target_atom: The atom to remove (optional).
    Returns:
        Nothing.
    """
    new_resi_name, target_proton = get_default_deproton_info(residue, target_atom)
    init_connectivity(target_atom, renew=False)
    if new_resi_name is None:
        if len(target_atom.attached_protons()) == 0:
            _LOGGER.info(f"target atom {target_atom} already have no H. keep original")
        else:
            _LOGGER.warning(f"cannot deprotonate {residue} on {target_atom}. keep original")
        return None
    _LOGGER.info(f"deprotonate {target_proton} from {residue}")
    residue.name = new_resi_name
    residue.find_atom_name(target_proton).delete_from_parent()
    #TODO rename/complete atoms after this (e.g.: ARG:NH1 case, HID -> HIE case)


def get_default_deproton_info(residue: Residue, target_atom: Union[None, Atom] = None) -> Tuple:
    """Returns the default proton in the residue on the target_atom (if provided) to deprotonate.
    Default HIP target is set to resulting HIE
    Args:
        residue: the target residue of deprotonation
        target_atom: the atom of deprotonation
    Return:
        a tuple of (residue_name_after_deproton, atom_name_of_the_removing_proton)
    """
    r_name = residue.name
    # default target atom
    if target_atom is None:
        target_atom_name = list(chem.residue.DEPROTONATION_MAPPER[r_name].keys())[0]
    else:
        target_atom_name = target_atom.name

    depro_info = chem.residue.DEPROTONATION_MAPPER.get(r_name, None)
    if depro_info is None:
        _LOGGER.warning(f"no default protonation info for {r_name}. Consider make a standard for it")
        return None, None
    if r_name in ["HIE", "HID"]:
        _LOGGER.warning(f"deprotonation info for {residue} is actually a switching between HID/HIE")
    target_atom_depro_info = depro_info.get(target_atom_name, None)
    if target_atom_depro_info is None:
        _LOGGER.warning(
            f"no default protonation info for {target_atom_name} in {r_name}. Could be no proton on it. Consider make a standard for it if do"
        )
        return None, None
    return target_atom_depro_info


def remove_side_chain_mutating_atom(residue: Residue, substitute_residue_name: str, remove_cb: bool = False) -> None:
    """remove mutating atoms of the side chain of the residue. The atoms are determined
    by the 3-letter name of the {substitute_residue_name}. (e.g.: to GLY requires removing more
    atoms)
    make changes in place.
    Args:
        residue:
            the target residue of change
        substitute_residue_name: 
            the 3-letter name of the residue after the subtitution mutation
        remove_cb:
            whether want to remove also CB when it is not removed"""
    # atoms to keep
    non_mutate_atom_names = chem.residue.get_non_mutate_atom_names(substitute_residue_name)
    if "CB" in non_mutate_atom_names and remove_cb:
        non_mutate_atom_names.remove("CB")

    # remove the rest atoms
    residue.remove_atoms_not_in_list(non_mutate_atom_names)


# == checker ==
def check_res_topology_error(
    stru: Structure,
    residue_key: Tuple[str, int],
    check_radius: float = 5.0,
):
    """check {stru} for topology error. (check for only the {residue_key} residue)
    i.e.: rings in structure should not be circling on other bonds.
    An example of this error is in https://github.com/ChemBioHTP/EnzyHTP/issues/110
    Args:
        stru: the target stru
        residue_key: key of the target residue
        check_radius: the radius of checking for the residue. (Unit: Ang)"""
    # 1. get connectivity for the structure
    # stru.init_connect()
    # TODO(qz)(high_prior)
    # 2. convert all bond in each residue to ploylines. find rings in each residue.
    # 3. check for any bond (polyline) from the target residue thread through any ring
    # 4. check for any ring from the target residue is threaded by any bond
    pass


def closest_n_residues( 
    target:Residue,
    n:int,
    method:str='centerofmass',
    include_H:bool=False) -> List[Residue]:
    target_coord = target.geom_center
    distances = list()
    for res in target.root().residues:
        if res == target:
            continue
        
        distances.append((mh.get_distance(target_coord, res.geom_center), res))

    distances = sorted(distances, key=lambda pr: pr[0])
    result = list()
    for (_, res) in distances[:n]:
        result.append( res )        
    return result
