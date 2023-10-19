"""Operations related to connectivity. Functions in this module initiate, modified, etc.
connectivity of Structure()

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2023-10-17
"""
# pylint: disable=function-redefined
from plum import dispatch
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
from enzy_htp import config as eh_config

# == API ==
@dispatch
def init_connectivity(target: Atom, renew: bool = False,
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
            _connect_res(res, maa_fix, ligand_fix, solvent_fix, metal_fix)
    else:
        _LOGGER.error(f"{target} dont have a parent residue, and thus connectivity is ill-defined")
        raise AttributeError

@dispatch
def init_connectivity(target: Residue, renew: bool = False,
        ligand_fix: str = "antechamber",
        maa_fix: str = "antechamber",
        solvent_fix: str = "caa",
        metal_fix: str = "isolate",):
    """init connectivty for Residue including all NCAAs."""
    if (not renew) and target.is_connected():
        return
    _connect_res(target, maa_fix, ligand_fix, solvent_fix, metal_fix)

@dispatch
def init_connectivity(target: Chain, renew: bool = False,
        ligand_fix: str = "antechamber",
        maa_fix: str = "antechamber",
        solvent_fix: str = "caa",
        metal_fix: str = "isolate",):
    """init connectivty for chain."""
    if (not renew) and target.is_connected():
        return

    for res in target:
        _connect_res(res, maa_fix, ligand_fix, solvent_fix, metal_fix)

@dispatch
def init_connectivity(target: Structure, renew: bool = False,
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
        _connect_res(res, maa_fix, ligand_fix, solvent_fix, metal_fix)

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
    if parent_residue.name in chem.solvent.RD_SOLVENT_LIST:
        cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP[parent_residue.name][atom.name]
    elif parent_residue.is_canonical():
        r = parent_residue
        r1 = parent_residue.chain[0]
        rm1 = parent_residue.chain[-1]
        if r is r1:
            # N terminal
            cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP_NTERMINAL[parent_residue.name][atom.name]
        else:
            if r == rm1:
                # C terminal
                cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP_CTERMINAL[parent_residue.name][atom.name]
            else:
                cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP[parent_residue.name][atom.name]
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
        metal_fix: str,) -> None:

    if res.is_canonical():
        _connect_caa(res)
    elif res.is_modified():
        _connect_maa(res, maa_fix)
    elif res.is_ligand():
        _connect_ligand(res, ligand_fix)
    elif res.is_solvent():
        _connect_solvent(res, solvent_fix)
    elif res.is_metal():
        _connect_metal(res, metal_fix)

def _connect_caa(res: Residue) -> None:
    """connect atoms in a canonical residue"""
    for atom in res.atoms:
        _connect_caa_atom(atom)

MOL_DESC_METHODS = ["antechamber"]
"""the list for method keywords that belongs to molecule description file based methods."""

MOL_DESC_GEN_MAPPER = {
    "antechamber" : None, # some func in amber interface
}
"""the mapper that interpret the method keyword to a function that generate mol describing file
The functions needs to support both ligand and modified amino acid"""

MOL_DESC_PARSER_MAPPER = {
    ".mol2" : Mol2Parser.get_structure,
    ".prepin" : PrepinParser.get_structure,
    ".prepi" : PrepinParser.get_structure,
}
"""the mapper that interpret the method keyword to a function that generate mol describing file"""

def _connect_ligand(lig: Ligand, method: str) -> None:
    """initiate connectivity for ligand"""
    support_method_list = ["antechamber"]
    if method in MOL_DESC_METHODS:
        _mol_desc_based_ncaa_method(lig, method)
        return

    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

def _connect_maa(maa: ModifiedResidue, method: str) -> None:
    """initiate connectivity for modified residue"""
    support_method_list = ["antechamber"]
    if method in MOL_DESC_METHODS:
        _mol_desc_based_ncaa_method(maa, method)
        return

    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

def _mol_desc_based_ncaa_method(ncaa: NonCanonicalBase, engine: str):
    """mol desc based method to generate connectivty for ncaa"""
    # 0. search lib for prepin/mol2 of maa
    ncaa_lib = eh_config["system.NCAA_LIB_PATH"]
    fs.safe_mkdir(ncaa_lib)
    mol_desc_path = search_ncaa_parm_file(ncaa,
                        target_method="any",
                        ncaa_lib_path=ncaa_lib)[0]

    # 1. make mol describing file for maa
    if not mol_desc_path:
        mol_desc_path = MOL_DESC_GEN_MAPPER[engine](ncaa)
    # 2. parse mol describing and clone into connectivity
    cnt_stru = MOL_DESC_PARSER_MAPPER[fs.get_file_ext(mol_desc_path)](mol_desc_path)
    cnt_maa = cnt_stru.residues[0]
    ncaa.clone_connectivity(cnt_maa)

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

def _connect_metal(met: MetalUnit, method: str) -> None:
    """initiate connectivity for metal"""
    support_method_list = ["isolate"]
    if method == "isolate":
        met.atom.connect = []

    if method not in support_method_list:
        _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

# def init_connect_for_ligands(self, method: str):
#     """"""
#     # TODO(qz)(high_prior)
#     # TODO generate prepi by itself and store it to a global database path so that other
#     # process in the same workflow can share the generated file.
#     support_method_list = ["caa"]
#     # ligand
#     # init
#     for lig in self.ligands:
#         for atom in lig:
#             atom.connect = []
#     # fix 1
#     if ligand_fix == 1:
#         for lig in self.ligands:
#             # read prepin for each ligand
#             with open(prepi_path[lig.name]) as f:
#                 line_id = 0
#                 if_loop = 0
#                 for line in f:
#                     line_id += 1
#                     if line.strip() == "":
#                         if if_loop == 1:
#                             # switch off loop and break if first blank after LOOP encountered
#                             if_loop = 0
#                             break
#                         continue
#                     if if_loop:
#                         lp = line.strip().split()
#                         lig._find_atom_name(lp[0]).connect.append(lig._find_atom_name(lp[1]))
#                         continue
#                     # loop connect starts at LOOP
#                     if line.strip() == "LOOP":
#                         if_loop = 1
#                         continue
#                     # coord starts at 11th
#                     if line_id >= 11:
#                         lp = line.strip().split()
#                         atom_id = str(lp[0]) - 3
#                         atom_cnt = str(lp[4]) - 3
#                         if atom_cnt != 0:
#                             lig[atom_id - 1].connect.append(lig[atom_cnt - 1])
#                             lig[atom_cnt - 1].connect.append(lig[atom_id - 1])
#     if method not in support_method_list:
#         _LOGGER.error(f"Method {method} not in supported list: {support_method_list}")

