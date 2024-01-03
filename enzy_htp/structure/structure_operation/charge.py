"""Operations related to charge. Functions in this module initiate, modified, etc.
atomic charges of Structure().

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2024-1-1
"""
# pylint: disable=function-redefined
from plum import dispatch
from typing import Union
import sys
from functools import partial

from ...structure import (
    Structure,
    Chain,
    Residue,
    NonCanonicalBase,
    Atom,
    Ligand,
    ModifiedResidue,
    MetalUnit,
    Solvent,
    Mol2Parser,
    PrepinParser,
    StructureRegion,
)

import enzy_htp.core.file_system as fs
from enzy_htp.core.logger import _LOGGER
from enzy_htp.core.exception import ResidueDontHaveAtom
import enzy_htp.chemical as chem
from enzy_htp._interface.ncaa_library import search_ncaa_parm_file
from enzy_htp import interface
from enzy_htp import config as eh_config

# == API ==
@dispatch
def init_charge(target: Atom, renew: bool = False, ncaa_lib: Union[str, None] = None,
        caa_fix: str = "ff19sb",
        ligand_fix: str = "bcc",
        maa_fix: str = "bcc",
        solvent_fix: str = "ff19sb",
        metal_fix: str = "formal_charge",
        ):
    """init charge for atom.
    Slow and redundant if applied to too many Atom().
    Use Structure or Residue as input instead in those cases"""
    if (not renew) and target.has_init_charge():
        return

    res: Residue = target.residue

    if res:
        if res.is_canonical():
            _init_charge_caa_atom(target)
        else:
            _init_charge_res(res, caa_fix, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)
    else:
        _LOGGER.error(f"{target} dont have a parent residue, and thus charge is ill-defined")
        raise AttributeError

@dispatch
def init_charge(target: Residue, renew: bool = False, ncaa_lib: Union[str, None] = None,
        caa_fix: str = "ff19sb",
        ligand_fix: str = "bcc",
        maa_fix: str = "bcc",
        solvent_fix: str = "ff19sb",
        metal_fix: str = "formal_charge",
    ):
    """init charge for Residue including all NCAAs."""
    if (not renew) and target.has_init_charge():
        return
    _init_charge_res(target, caa_fix, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)

@dispatch
def init_charge(target: Chain, renew: bool = False, ncaa_lib: Union[str, None] = None,
        caa_fix: str = "ff19sb",
        ligand_fix: str = "bcc",
        maa_fix: str = "bcc",
        solvent_fix: str = "ff19sb",
        metal_fix: str = "formal_charge",
    ):
    """init charge for chain."""
    if (not renew) and target.has_init_charge():
        return

    for res in target:
        _init_charge_res(res,caa_fix, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)

@dispatch
def init_charge(target: StructureRegion, renew: bool = False, ncaa_lib: Union[str, None] = None,
        caa_fix: str = "ff19sb",
        ligand_fix: str = "bcc",
        maa_fix: str = "bcc",
        solvent_fix: str = "ff19sb",
        metal_fix: str = "formal_charge",
        cap_fix: str = "fixed_charge",
    ) -> None:
    for res in target.involved_residues:
        if (not renew) and res.has_init_charge():
            continue
        _init_charge_res(res, caa_fix, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)

    if cap_fix == "fixed_charge":
        for cap in target.caps:
            cap.init_fixed_atomic_charges()
    else:
        _LOGGER.error(f"cap_fix = {cap_fix} not supported.")
        raise ValueError

    # check if all atoms have charge
    for atm in target.atoms:
        if not atm.has_init_charge():
            _LOGGER.error(f"Atom {atm} doesn't have charge record after initiation.")
            sys.exit(1)

@dispatch
def init_charge(target: Structure, renew: bool = False, ncaa_lib: Union[str, None] = None,
        caa_fix: str = "ff19sb",
        ligand_fix: str = "bcc",
        maa_fix: str = "bcc",
        solvent_fix: str = "ff19sb",
        metal_fix: str = "formal_charge",
    ) -> None:
    """
    Initiate charge for the Structure. (dispatch available for Chain, Residue, Atom)
    Save the charge to self._connect of each Atom().
    **  This function is placed here instead of inside Structure to use _strerface for 
        charge generation of NCAAs

    Args:
        renew:
            whether renew the existing charge
        ncaa_lib:
            the path of ncaa library for ligand/maa/metal_fix. (default: eh_config["system.NCAA_LIB_PATH"])
        ligand_fix:
            the method that determines charge for ligand. (see details below)
        metal_fix:
            the method that determines charge for metal. (see details below)
        maa_fix:
            the method that determines charge for modified animo acid. (see details below)
        solvent_fix:
            the method that determines charge for solvent. (see details below)

    Details (method for connecting each component of each keyword):
        polypeptide:
            ff19sb:
                using documented charge for each canonical amino acid from ff19sb.
        ligand:
            fix = "bcc":
                use bcc charge
        modified residues:
            fix = "bcc": same as above in ligand.
        solvent:
            fix = "ff19sb": same as polypeptide part
            TODO make new fix for non-water solvent when encountered
        metalatom:
            fix = "formal_charge": just use formal charge from user define.
            fix = "mcpb": connect to donor atom (MCPB?) TODO
    """
    for res in target.residues:
        if (not renew) and res.has_init_charge():
            continue
        _init_charge_res(res,caa_fix, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)

    # check if all atoms have charge
    for atm in target.atoms:
        if not atm.has_init_charge():
            _LOGGER.error(f"Atom {atm} doesn't have charge record after initiation.")
            sys.exit(1)

# == parts ==
def _init_charge_caa_atom(atom: Atom, force_field: str) -> None:
    """
    Initiate charge for {atom} in a canonical amino acid and water
    find connect atom base on:
    1. chem.residue.RESIDUE_CHARGE_MAP[force_field]
    2. parent residue name
    * Using standard Amber atom names and C/N terminal name.
        (TODO make this a standard and convert other atom name formats)
    save found charge value to atom._charge
    """
    # san check
    if force_field not in chem.residue.RESIDUE_CHARGE_MAP:
        _LOGGER.error(f"given force_field ({force_field}) not supported in this function. Supported: {chem.residue.RESIDUE_CHARGE_MAP.keys()}")

    parent_residue = atom.parent
    res_name = parent_residue.name
    if res_name == "HIS": # deal with HIS
        _LOGGER.warning("HIS found in the target structure, treat as HIE by default."
                        "(consider protonate the structure first using enzy_htp.preparation.protonate_stru())")
        res_name = "HIE"
    if res_name in chem.solvent.RD_SOLVENT_LIST:
        atom_charge = chem.residue.RESIDUE_CHARGE_MAP[force_field][res_name][atom.name]
    elif parent_residue.is_canonical():
        r = parent_residue
        r1 = parent_residue.chain.n_ter_residue()
        rm1 = parent_residue.chain.c_ter_residue()
        if r is r1:
            # N terminal
            atom_charge = chem.residue.RESIDUE_CHARGE_MAP_NTERMINAL[force_field][res_name][atom.name]
        else:
            if r == rm1:
                # C terminal
                atom_charge = chem.residue.RESIDUE_CHARGE_MAP_CTERMINAL[force_field][res_name][atom.name]
            else:
                atom_charge = chem.residue.RESIDUE_CHARGE_MAP[force_field][res_name][atom.name]
    else:
        _LOGGER.error(
            f"wrong method of getting charge of non-canonical residue {atom.parent}. Use _init_charge_maa or _init_charge_ligand etc.")
        sys.exit(1)

    atom.charge = atom_charge

def _init_charge_res(
        res: Residue,
        caa_fix: str,
        maa_fix: str,
        ligand_fix: str,
        solvent_fix: str,
        metal_fix: str,
        ncaa_lib: str) -> None:

    if res.is_canonical():
        _init_charge_caa(res, caa_fix)
    elif res.is_modified():
        _init_charge_maa(res, maa_fix, ncaa_lib)
    elif res.is_ligand():
        _init_charge_ligand(res, ligand_fix, ncaa_lib)
    elif res.is_solvent():
        _init_charge_solvent(res, solvent_fix)
    elif res.is_metal():
        _init_charge_metal(res, metal_fix, ncaa_lib)

def _init_charge_caa(res: Residue, method: str) -> None:
    """connect atoms in a canonical residue"""
    support_method_list = ["ff19sb"]
    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

    if method == "ff19sb":
        for atom in res.atoms:
            _init_charge_caa_atom(atom, method)
    else:
        # TODO none fixed charge method
        pass        

CHG_METHOD_EXT_MAPPER = {
    "bcc" : "prepin",
    "resp" : "gout",
}
"""the mapper that maps charge methods names to the extension used in the engine"""

CHG_ENGINE_MAPPER = {
    "bcc" : partial(interface.amber.antechamber_ncaa_to_moldesc, charge_method="AM1BCC"),
    "resp" : interface.gaussian.build_single_point_engine, # TODO make some partial
}
"""the mapper that interpret the method keyword to a function that generate mol charge file
The functions needs to support both ligand and modified amino acid"""

CHG_METHODS = ["bcc", "resp"]
"""the list for method keywords that indicates a charge model of support."""

CHG_FILE_PARSER_MAPPER = {
    ".mol2" : Mol2Parser.get_structure,
    ".prepin" : PrepinParser.get_structure,
    ".gout" : None, # TODO make this
}
"""the mapper that finds the parser of a specific file type"""

def _init_charge_ligand(lig: Ligand, method: str, ncaa_lib: str) -> None:
    """initiate charge for ligand"""
    support_method_list = ["bcc", "resp"]
    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

    _mol_desc_based_ncaa_method(lig, method, ncaa_lib)

def _init_charge_maa(maa: ModifiedResidue, method: str, ncaa_lib: str) -> None:
    """initiate charge for modified residue"""
    raise Exception("TODO") # TODO look into whether .ac is just enough
    support_method_list = ["antechamber"]
    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")
    _mol_desc_based_ncaa_method(maa, method, ncaa_lib)

def _mol_desc_based_ncaa_method(ncaa: NonCanonicalBase, method: str, ncaa_lib: str):
    """mol desc based method to generate charge for ncaa"""
    # 0. search lib for charge files of ncaa
    if ncaa_lib is None:
        ncaa_lib = eh_config["system.NCAA_LIB_PATH"]
    fs.safe_mkdir(ncaa_lib)
    mol_desc_path = search_ncaa_parm_file(
                        ncaa,
                        target_method=f"CHGONLY-{method.upper()}",
                        ncaa_lib_path=ncaa_lib)[0]

    # 1. make mol describing file for ncaa
    if not mol_desc_path:
        ext = CHG_METHOD_EXT_MAPPER[method]
        mol_desc_path = f"{ncaa_lib}/{ncaa.name}_CHGONLY-{method.upper()}.{ext}"
        CHG_ENGINE_MAPPER[method](ncaa=ncaa, out_path=mol_desc_path)

    # 2. parse mol describing and clone into charge
    cnt_stru = CHG_FILE_PARSER_MAPPER[fs.get_file_ext(mol_desc_path)](mol_desc_path)
    cnt_ncaa = cnt_stru.residues[0]
    ncaa.clone_charge(cnt_ncaa) # TODO make this

def _init_charge_solvent(sol: Solvent, method: str) -> None:
    """initate charge for solvent."""
    support_method_list = ["ff19sb"]
    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")
        raise ValueError

    if method == "ff19sb":
        for atom in sol.atoms:
            _init_charge_caa_atom(atom, method)

def _init_charge_metal(met: MetalUnit, method: str, ncaa_lib: str) -> None:
    """initiate charge for metal
    formal_charge:
        just use formal charge of the metal from user assignment"""
    raise Exception("TODO")
    support_method_list = ["formal_charge"]
    if method == "isolate":
        met.atom.connect = []

    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")
