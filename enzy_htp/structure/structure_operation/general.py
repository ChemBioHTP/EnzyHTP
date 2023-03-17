"""General opertion of Structure(). functions in this module take Structure() as input and do 
operations on it. Place holder for uncategorized functions.

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2022-09-19
"""

import copy
from typing import Tuple, Union
from plum import dispatch

from enzy_htp.core.logger import _LOGGER
import enzy_htp.chemical as chem
from ..structure import Structure, Solvent, Chain, Residue, Atom


def remove_solvent(stru: Structure) -> Structure:
    """
    remove all Solvent() for {stru}.
    Make changes in-place and return a reference of the changed
    original object.
    """
    _LOGGER.debug(f"removing {len(stru.solvents)} solvents")
    solv: Solvent
    for solv in stru.solvents:
        solv.delete_from_parent()

    return stru


def remove_empty_chain(stru: Structure) -> Structure:
    """
    remove empty chains
    Make changes in-place and return a reference of the changed
    original object.
    """
    ch: Chain
    for ch in copy.copy(stru.chains):
        if ch.is_empty():
            _LOGGER.debug(f"removing {ch}")
            ch.delete_from_parent()
    return stru


def remove_non_peptide(stru: Structure) -> Structure:
    """remove the non-peptide parts of the structure. 
    Make changes in-place and return a reference of the changed original object."""
    non_peptides = list(filter(lambda c: not c.is_polypeptide(), stru.chains))
    ch: Chain
    for ch in non_peptides:
        ch.delete_from_parent()
    return stru


@dispatch
def update_residues(stru: Structure, ref_stru: Structure) -> Structure:
    """
    Update atoms and residue names to residues in the stru
    The sequence should holds constant since it serves as reference
    Args:
        stru: the target structure
        ref_stru: the reference structure
    Returns:
        stru: the changed original structure
    """
    #san check for ref_stru
    stru.is_idx_subset(ref_stru)

    # update residues in stru with correponding residue idxes in ref_stru
    self_res_mapper = stru.residue_mapper
    ref_res: Residue
    for ref_res in ref_stru.residues:
        self_res = self_res_mapper[ref_res.key()]
        if self_res.name != ref_res.name:
            _LOGGER.info(f"updating {self_res.key()} {self_res.name} to {ref_res.name}")
            self_res.name = ref_res.name
        self_res.atoms = copy.deepcopy(ref_res.atoms)  # this will also set self_res as parent
    stru.renumber_atoms()
    return stru


@dispatch
def update_residues(resi: Residue, ref_resi: Residue) -> Residue:  # pylint: disable=function-redefined
    """
    Update atoms and residue names to the single residue
    (it doesnt matter if there are same in sequence)
    Args:
        resi: the target residue
        ref_resi: the reference residue
    Return:
        the changed original residue
    """
    if resi.name != ref_resi.name:
        _LOGGER.info(f"updating {resi.key()} {resi.name} to {ref_resi.name}")
        resi.name = ref_resi.name
    resi.atoms = copy.deepcopy(ref_resi.atoms)  # this will also set resi as parent
    return resi


def deprotonate_residue(residue: Residue, target_atom: Union[None, Atom] = None) -> None:
    """
    deprotonate the {residue} on the {target_atom} (if provided).
    remove the acidic hydrogen attached to the {target_atom} and change the residue name
    correponding to chem.residue.DEPROTONATION_MAPPER or /resource/ProtonationState.cdx
    """
    new_resi_name, target_proton = get_default_deproton_info(residue, target_atom)
    if new_resi_name is None:
        if len(target_atom.attached_protons()) == 0:
            _LOGGER.info(f"target atom {target_atom} already have no H. keep original")
        else:
            _LOGGER.warning(
                f"cannot deprotonate {residue} on {target_atom}. keep original")
        return None
    _LOGGER.info(f"deprotonate {target_proton} from {residue}")
    residue.name = new_resi_name
    residue.find_atom_name(target_proton).delete_from_parent()
    #TODO rename/complete atoms after this (e.g.: ARG:NH1 case, HID -> HIE case)


def get_default_deproton_info(residue: Residue,
                              target_atom: Union[None, Atom] = None) -> Tuple:
    """
    return the default proton in the residue on the target_atom (if provided) to deprotonate
    Default HIP target is set to resulting HIE
    """
    r_name = residue.name
    # default target atom
    if target_atom is None:
        target_atom_name = list(chem.residue.DEPROTONATION_MAPPER[r_name].keys())[0]
    else:
        target_atom_name = target_atom.name

    depro_info = chem.residue.DEPROTONATION_MAPPER.get(r_name, None)
    if depro_info is None:
        _LOGGER.warn(
            f"no default protonation info for {r_name}. Consider make a standard for it")
        return None, None
    if r_name in ["HIE", "HID"]:
        _LOGGER.warn(
            f"deprotonation info for {residue} is actually a switching between HID/HIE")
    target_atom_depro_info = depro_info.get(target_atom_name, None)
    if target_atom_depro_info is None:
        _LOGGER.warn(
            f"no default protonation info for {target_atom_name} in {r_name}. Could be no proton on it. Consider make a standard for it if do"
        )
        return None, None
    return target_atom_depro_info
