"""Operations related to connectivity. Functions in this module initiate, modified, etc.
connectivity of Structure()

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-10-17
"""
# pylint: disable=function-redefined
from enzy_htp.structure.structure_region.api import create_region_from_residues
from enzy_htp.structure.structure_region.structure_region import StructureRegion
from plum import dispatch
from typing import Union
import sys

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
def init_connectivity(target: Atom, renew: bool = False, ncaa_lib: Union[str, None] = None,
        ligand_fix: str = "antechamber",
        maa_fix: str = "antechamber",
        solvent_fix: str = "caa",
        metal_fix: str = "isolate",
        ):
    """init connectivty for atom.
    Slow and redundant if applied to too many Atom().
    Use Structure or Residue as input instead in those cases"""
    if (not renew) and target.is_connected():
        return

    res: Residue = target.residue

    if res:
        if res.is_canonical():
            _connect_caa_atom(target)
        else:
            _connect_res(res, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)
    else:
        _LOGGER.error(f"{target} dont have a parent residue, and thus connectivity is ill-defined")
        raise AttributeError

@dispatch
def init_connectivity(target: Residue, renew: bool = False, ncaa_lib: Union[str, None] = None,
        ligand_fix: str = "antechamber",
        maa_fix: str = "antechamber",
        solvent_fix: str = "caa",
        metal_fix: str = "isolate",):
    """init connectivty for Residue including all NCAAs."""
    if (not renew) and target.is_connected():
        return
    _connect_res(target, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)

@dispatch
def init_connectivity(target: Chain, renew: bool = False, ncaa_lib: Union[str, None] = None,
        ligand_fix: str = "antechamber",
        maa_fix: str = "antechamber",
        solvent_fix: str = "caa",
        metal_fix: str = "isolate",):
    """init connectivty for chain."""
    if (not renew) and target.is_connected():
        return

    for res in target:
        _connect_res(res, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)

@dispatch
def init_connectivity(target: Structure, renew: bool = False, ncaa_lib: Union[str, None] = None,
        ligand_fix: str = "antechamber",
        maa_fix: str = "antechamber",
        solvent_fix: str = "caa",
        metal_fix: str = "isolate") -> None:
    """
    Initiate connectivity for the Structure. (dispatch available for Chain, Residue, Atom)
    Save the connectivity to self._connect of each Atom().
    **  This function is placed here instead of inside Structure to use _strerface for 
        connectivity generation of NCAAs

    Args:
        renew:
            whether renew the existing connectivity
        ncaa_lib:
            the path of ncaa library for ligand/maa/metal_fix. (default: eh_config["system.NCAA_LIB_PATH"])
        ligand_fix:
            the method that determines connectivity for ligand. (see details below)
        metal_fix:
            the method that determines connectivity for metal. (see details below)
        maa_fix:
            the method that determines connectivity for modified animo acid. (see details below)
        solvent_fix:
            the method that determines connectivity for solvent. (see details below)

    Details (method for connecting each component of each keyword):
        polypeptide:
            using documented connectivity for each canonical amino acid from Amber library.
        ligand:
            fix = "antechamber":
                use antechamber to generated connectivity and read from prepin file.
                (according to https://ambermd.org/doc/prep.html the coordniate line will
                always start at the 11th line after 3 DUMM.)
        modified residues:
            fix = "antechamber": same as above in ligand.
        solvent:
            fix = "caa": same as polypeptide part
            TODO make new fix for non-water solvent when encountered
        metalatom:
            fix = "isolate": treat as isolated atom
            fix = "mcpb": connect to donor atom (MCPB?) TODO
    """
    for res in target.residues:
        if (not renew) and res.is_connected():
            continue
        _connect_res(res, maa_fix, ligand_fix, solvent_fix, metal_fix, ncaa_lib)

    # check if all atoms are connected
    for atm in target.atoms:
        if not atm.is_connected():
            _LOGGER.error(f"Atom {atm} doesn't have connect record after initiation.")
            sys.exit(1)

# == parts ==
def _connect_caa_atom(atom: Atom) -> None:
    """
    Initiate connectivity for {atom} in a canonical amino acid and water
    find connect atom base on:
    1. chem.residue.RESIDUE_CONNECTIVITY_MAP
    2. parent residue name
    * Using standard Amber atom names and C/N terminal name.
        (TODO make this a standard and convert other atom name formats)
    save found list of Atom object to atom._connect
    """
    connect = []
    parent_residue = atom.parent
    res_name = parent_residue.name
    if res_name == "HIS": # deal with HIS
        _LOGGER.warning("HIS found in the target structure, treat as HIE by default."
                        "(consider protonate the structure first using enzy_htp.preparation.protonate_stru())")
        res_name = "HIE"
    if res_name in chem.solvent.RD_SOLVENT_LIST:
        cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP[res_name][atom.name]
    elif parent_residue.is_canonical():
        r = parent_residue
        r1 = parent_residue.chain.n_ter_residue()
        rm1 = parent_residue.chain.c_ter_residue()
        if r is r1:
            # N terminal
            cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP_NTERMINAL[res_name][atom.name]
        else:
            if r == rm1:
                # C terminal
                cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP_CTERMINAL[res_name][atom.name]
            else:
                cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP[res_name][atom.name]
    else:
        _LOGGER.error(
            f"wrong method of getting connectivity of non-canonical residue {atom.parent}. Use connect_maa or connect_ligand etc.")
        sys.exit(1)

    for name in cnt_atomnames:
        try:
            if name not in ["-1C", "+1N"]:
                cnt_atom = parent_residue.find_atom_name(name)
            if name == "-1C":
                cnt_resi = parent_residue.chain.find_residue_idx(parent_residue.idx - 1)
                cnt_atom = cnt_resi.find_atom_name("C")
            if name == "+1N":
                cnt_resi = parent_residue.chain.find_residue_idx(parent_residue.idx + 1)
                cnt_atom = cnt_resi.find_atom_name("N")
            connect.append((cnt_atom, None))
        except ResidueDontHaveAtom as e:
            _LOGGER.warning(f"missing connecting atom {e.atom_name} of {atom}. Structure maybe incomplete.")
    atom.connect = connect

def _connect_res(res: Residue,
        maa_fix: str,
        ligand_fix: str,
        solvent_fix: str,
        metal_fix: str,
        ncaa_lib: str) -> None:

    if res.is_canonical():
        _connect_caa(res)
    elif res.is_modified():
        _connect_maa(res, maa_fix, ncaa_lib)
    elif res.is_ligand():
        _connect_ligand(res, ligand_fix, ncaa_lib)
    elif res.is_solvent():
        _connect_solvent(res, solvent_fix)
    elif res.is_metal():
        _connect_metal(res, metal_fix, ncaa_lib)

def _connect_caa(res: Residue) -> None:
    """connect atoms in a canonical residue"""
    for atom in res.atoms:
        _connect_caa_atom(atom)

MOL_DESC_METHODS = ["antechamber"]
"""the list for method keywords that belongs to molecule description file based methods."""

MOL_DESC_GEN_MAPPER = {
    "antechamber" : interface.amber.antechamber_ncaa_to_moldesc,
}
"""the mapper that interpret the method keyword to a function that generate mol describing file
The functions needs to support both ligand and modified amino acid"""

MOL_DESC_PARSER_MAPPER = {
    ".mol2" : Mol2Parser.get_structure,
    ".prepin" : PrepinParser.get_structure,
    ".prepi" : PrepinParser.get_structure,
}
"""the mapper that interpret the method keyword to a function that generate mol describing file"""

def _connect_ligand(lig: Ligand, method: str, ncaa_lib: str) -> None:
    """initiate connectivity for ligand"""
    support_method_list = ["antechamber"]
    if method in MOL_DESC_METHODS:
        _mol_desc_based_ncaa_method(lig, method, ncaa_lib)
        return

    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

def _connect_maa(maa: ModifiedResidue, method: str, ncaa_lib: str) -> None:
    """initiate connectivity for modified residue"""
    support_method_list = ["antechamber"]
    if method in MOL_DESC_METHODS:
        _mol_desc_based_ncaa_method(maa, method, ncaa_lib)
        return

    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

def _mol_desc_based_ncaa_method(ncaa: NonCanonicalBase, engine: str, ncaa_lib: str):
    """mol desc based method to generate connectivty for ncaa"""
    # 0. search lib for prepin/mol2 of ncaa
    if ncaa_lib is None:
        ncaa_lib = eh_config["system.NCAA_LIB_PATH"]
    fs.safe_mkdir(ncaa_lib)
    mol_desc_path = search_ncaa_parm_file(ncaa,
                        target_method="any",
                        ncaa_lib_path=ncaa_lib)[0]

    # 1. make mol describing file for ncaa
    if not mol_desc_path:
        mol_desc_path = f"{ncaa_lib}/{ncaa.name}_any.prepin"  # swicth to mol2 after finish all unit tests of mol2_io
        MOL_DESC_GEN_MAPPER[engine](ncaa=ncaa, out_path=mol_desc_path)
    # 2. parse mol describing and clone into connectivity
    cnt_stru = MOL_DESC_PARSER_MAPPER[fs.get_file_ext(mol_desc_path)](mol_desc_path)
    cnt_ncaa = cnt_stru.residues[0]

    ncaa.clone_connectivity(cnt_ncaa)

def _connect_solvent(sol: Solvent, method: str) -> None:
    """initate connectivity for solvent."""
    support_method_list = ["caa"]
    if method == "caa":
        if sol.name not in chem.solvent.RD_SOLVENT_LIST:
            _LOGGER.error(f"'caa' method: solvent {sol.name} not in supported list: {chem.solvent.RD_SOLVENT_LIST}")
            raise ValueError
        for atom in sol.atoms:
            _connect_caa_atom(atom)

    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")
        raise ValueError

def _connect_metal(met: MetalUnit, method: str, ncaa_lib: str) -> None:
    """initiate connectivity for metal"""
    support_method_list = ["isolate"]
    if method == "isolate":
        met.atom.connect = []

    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")
