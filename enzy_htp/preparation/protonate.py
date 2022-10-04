"""Defines functions for protonating Structure(), Ligand() and MetalUnit() objects as well as raw PDB files.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-05
"""
# TODO(CJ): add more documentation
from pathlib import Path
from typing import Union,Tuple, Dict

import enzy_htp.core as core
from enzy_htp.core import file_system as fs
from enzy_htp import config
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import (
    Structure,
    Ligand,
    structure_from_pdb,
    ligand_from_pdb,
    Chain,
    PDBParser
)
from enzy_htp.structure.metal_atom import MetalUnit
import enzy_htp.structure.structure_operation as stru_oper

from pdb2pqr.main import main_driver as run_pdb2pqr
from pdb2pqr.main import build_main_parser as build_pdb2pqr_parser
import openbabel
import openbabel.pybel as pybel
from .pdb_line import read_pdb_lines

def protonate_stru(stru: Structure,
                   ph: float = 7.0, protonate_ligand: bool = False,
                   engine: str = "pdb2pqr", ligand_engine: str = "pybel",
                   **kwargs) -> Structure:
    """
    This science API solves the protein protonation problem.
    that add missing H atoms to the {stru}. Protonation states are determined 
    for residues with multiple ones.
    Args:
        stru: the input structure
        ph: the pH value for determining the protonation state
        protonate_ligand: if also protonate ligand
        engine: engine for determining the pKa and adding hydrogens to the protein peptide part
            (current available keywords):
            pdb2pqr
        ligand_engine: engine for adding hydrogens to ligands
            (current available keywords):
            pybel
        **kwarg: setting/option related to specific engine. TODO figure out a better way to doc this
    Returns:
        the reference of the in-place changed {stru}

    Details:
        Approximately 88% of the structures in the protein data bank (PDB) are determined
    by X-ray crystallography, which can not, in general, resolve positions of most hydrogen
    atoms. The same problem appear in structures obtained from structure prediction tools
    too, AlphaFold2 also cannot give accurate position for hydrogens.
    (https://github.com/deepmind/alphafold/issues/598) Thus accuately determine protonation
    state is a vital part of structural preparation in EnzyHTP to ensure the accuracy of the
    modeling result.
    In short, the challenge is to predicting the protonation states of titratable groups such
    as the side chains of ASP, GLU, ARG, LYS, TYR, HIS, CYS, and ligands.

    Avaible strageties of protonation in the field are: (TODO do more search on this)
    Protein Protonation:
        Emperical pKa Prediction:
        - PROPKA3 (https://pubs.acs.org/doi/10.1021/ct100578z, https://github.com/jensengroup/propka)
            PROPKA3 trained parameters for a complex energy function. Commonly used in the field.
            Consider H-bonding network optimization.
            Consider considering H-bond with ligand.
            Not consider coordination to metal (https://github.com/jensengroup/propka/issues/135)
        - Rosetta-pH (https://www.sciencedirect.com/science/article/pii/S0006349512007333)
            a augmented scoring function from REF was created for pKa prediction.

        Poisson Boltzman (PB) model
        - H++ (https://academic.oup.com/nar/article/40/W1/W537/1072301, http://newbiophysics.cs.vt.edu/H++/)
            The only software that can consider the coordination to metal. Most commonly used in the field
            Consider coordination to metal.
            Unknow TODO H-bonding network optimization.
            Unknow TODO considering H-bond with ligand.
        - PypKa (https://pubs.acs.org/doi/full/10.1021/acs.jcim.0c00718, https://github.com/mms-fcul/PypKa)
            TODO add a summary
            Not consider coordination to metal (https://github.com/mms-fcul/PypKa/issues/6)

        Constant pH MD (believe to be most accurate but considering slow)
        - Amber (http://ambermd.org/tutorials/advanced/tutorial18/index.htm)

        There is also a experimently determined database: 
        - PKAD https://academic.oup.com/database/article/doi/10.1093/database/baz024/5359213
    Ligand Protonation:
        TODO need to figure out the algrothim #12
        - OpenBable/PyBel (Default)
            often have poor accuracy. 
        - reduce (http://kinemage.biochem.duke.edu/software/README.reduce.html) TODO
        - Dimorphite (from https://durrantlab.pitt.edu/dimorphite-dl/) TODO
        - Protons from OpenMM (https://protons.readthedocs.io/en/latest/) TODO
        - MCCE2 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2735604/) TODO
        - Protonate3D (http://www.ccl.net/cca/documents/proton/proton.htm) TODO
    """
    PEPTIDE_PROTONATION_METHODS[engine](stru, ph, **kwargs)
    if protonate_ligand:
        LIGAND_PROTONATION_METHODS[ligand_engine](stru, ph, **kwargs)

def protonate_peptide_with_pdb2pqr(stru: Structure, ph: float = 7.0,
                                   int_pdb_path: Union[str, None] = None,
                                   int_pqr_path: Union[str, None] = None,
                                   metal_fix_method: str = "deprotonate_all"):
    """
    Add missing hydrogens and determine protonation state of the peptide part of protein
    using [PDB2PQR](https://www.poissonboltzmann.org/) via the pdb2pqr python [package](https://pdb2pqr.readthedocs.io/en/latest/).
    Change the stru in place. (alignment of PDB2PQR result and original stru is based on residue keys)
    TODO(qz): support in-python implemtation of this
    Args:
        stru: the target Structure()
        ph: the target pH
        int_pqr_path: path for intermediate pqr file (not this will be changed to pdb extension)
        int_pdb_path: path for intermediate pdb file
    Returns:
        stru: a reference of the changed original structure
    """
    sp = PDBParser()
    if int_pdb_path is None:
        fs.safe_mkdir(config["system.SCRATCH_DIR"]) # make them together into make_temp_file
        int_pdb_path = fs.get_valid_temp_name(f"{config['system.SCRATCH_DIR']}/protonate_peptide_with_pdb2pqr_input.pdb")
    if int_pqr_path is None:
        fs.safe_mkdir("./temp")
        int_pqr_path = fs.get_valid_temp_name(f"{config['system.SCRATCH_DIR']}/protonate_peptide_with_pdb2pqr_output.pdb")
    if fs.get_file_ext(int_pqr_path) == ".pqr":
        _LOGGER.warning(f"changing {int_pqr_path} extension to pdb. This filename now changes.")
        int_pqr_path = fs.get_valid_temp_name(int_pqr_path.removesuffix("pqr")+"pdb")
        if int_pqr_path == int_pdb_path:
            _LOGGER.warning("int_pqr_path and int_pdb_path is the same after extension change. Adding an index.")
            int_pqr_path = fs.get_valid_temp_name(int_pqr_path.removesuffix(".pdb")+"_1.pdb")

    with open(int_pdb_path, "w") as of:
        of.write(sp.get_file_str(stru)) # give the whole structure as input here as PropKa can use ligand
    pdb2pqr_protonate_pdb(int_pdb_path, int_pqr_path, ph)
    peptide_protonated_stru = sp.get_structure(int_pqr_path)
    stru_oper.remove_non_peptide(peptide_protonated_stru) # keep the peptide only (sometime it has solvent)
    stru_oper.update_residues(stru, peptide_protonated_stru)
    protonate_peptide_fix_metal_donor(stru, method=metal_fix_method)
    if _LOGGER.level > 10: # not DEBUG or below
        fs.safe_rm(int_pdb_path)
        fs.safe_rm(int_pqr_path)
    return stru

def pdb2pqr_protonate_pdb(
    pdb_path: str, pqr_path: str, ph: float = 7.0, ffout: str = "AMBER"
) -> None:
    """
    This is warpper function of pdb2pqr.
    Runs PDB2PQR on a specified pdb file and saves it to the specified pqr path. This preparation step
    makes use of [PDB2PQR](https://www.poissonboltzmann.org/) via the pdb2pqr python [package](https://pdb2pqr.readthedocs.io/en/latest/).
    Adds in missing atoms and finds the protonation state of the pdb file.
    """
    # TODO(CJ): check if ffout is valid.
    # TODO(CJ): maybe improve the documentation here?
    core.check_valid_ph(ph)
    pdb2pqr_parser = build_pdb2pqr_parser()
    args = pdb2pqr_parser.parse_args(
        [
            "--ff=PARSE",
            "--ffout=" + ffout,
            "--with-ph=" + str(ph),
            "--log-level=CRITICAL",
            pdb_path,
            pqr_path,
        ]
    )
    core._LOGGER.info(f"Running pdb2pqr on '{pdb_path}'...")
    run_pdb2pqr(args)
    core._LOGGER.info(f"Finished running pdb2pqr! Output saved to '{pqr_path}'")

def protonate_peptide_fix_metal_donor(stru: Structure, method="deprotonate_all"):
    """
    fix protonation state around the metal coordination center
    Args:
        stru: target structure. Should be after peptide protonation first
        method: method for determine the protonation state of donor residues
            (current available keywords)
            deprotonate_all: 
    Return:
        (change stru in place)
        a reference of the changed stru
    """
    center: MetalUnit
    for center in stru.metalcenters:
        METAL_FIX_METHODS[method](center)

def deprotonate_metal_donors(center: MetalUnit):
    """
    deprotonate all donor atoms from donor residues away from neutral.
    apply change to the parent structure of this metal {center}
    """
    donor_mapper = center.get_donor_mapper(method = "ionic")
    for d_resi, d_atoms in donor_mapper.items():
        if d_resi.is_deprotonatable(): # the donor atom selection guarantees the atom is deprotonable
            # find_closest_h_to_center(d_atom, center)
            stru_oper.deprotonate_residue(d_resi, d_atoms[0]) # TODO(qz): refine this by also determine the closest proton
        elif d_resi.is_hetatom_noproton():
            _LOGGER.info(f"donor residue {d_resi} already have no proton in center {center}")
        else:
            _LOGGER.warn(f"uncommon donor residue {d_resi} found in center {center}")

def protonate_metal_donors_with_pka_recalc(center: MetalUnit):
    """place holder for most accurate metal fix method"""
    pass

METAL_FIX_METHODS = {
    "deprotonate_all": deprotonate_metal_donors
}

def protonate_ligand_with_pybel(stru: Structure, ph: float = 7.0, int_ligand_file_path = None):
    """
    the inteface for using PYBEL to protonate all ligands in {stru} with a given ph
    Args:
        stru: the target stru to protonate
        """

def pybel_protonate_pdb_ligand(in_path: str, out_path: str, ph: float = 7.0):
    """
    This is a wrapper of [PYBEL](https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html)
    that takes a ligand in PDB format and add missing H atoms with the given pH in the output file
    """
    outp1_path = path[:-4]+'_badname_aH.pdb'
    out_path = path[:-4]+'_aH.pdb'

    pybel.ob.obErrorLog.SetOutputLevel(0)
    mol = next(pybel.readfile('pdb', in_path))
    mol.OBMol.AddHydrogens(False, True, ph)
    mol.write('pdb', outp1_path, overwrite=True)
    # fix atom label abd determing net charge
    if keep_name:
        cls._fix_ob_output(outp1_path, out_path, ref_name_path=in_path)
    else:
        cls._fix_ob_output(outp1_path, out_path)
    return out_path



def protonate_missing_elements(
    old_pdb_path: str, new_pdb_path: str, work_dir: str
) -> Structure:
    """Method that compares two Structure() objects and combines missin elements from old_str to new_stru"""
    # TODO(CJ) Maybe this should be part of the structure.structure.py file
    old_stru = structure_from_pdb(old_pdb_path)
    new_stru = structure_from_pdb(new_pdb_path)

    if old_stru == new_stru:  # CHECK(CJ): Only sequence elements.
        return new_stru
        new_stru = self.merge_structure_elements(old_stru, new_stru)

    # protonate ligands and combine with the pqr file
    ligand_list = old_stru.ligands
    if len(ligand_list):
        # print('lig-list',ligand_list)
        core._LOGGER.info(f"Merging {len(ligand_list)} ligands in old structure!")
        lig_dir = work_dir + "/ligands/"
        fs.safe_mkdir(lig_dir)
        # TODO: logging
        # print(lig_dir)
        # print("-asdgasdgasg")
        old_keys = list(map(lambda l: l.residue_key, ligand_list))
        new_ligands = list(
            map(
                lambda ll: protonate_ligand(ll, dirname=lig_dir, ph=7),
                ligand_list,
            )
        )
        # print(new_ligands)
        for lig, ok in zip(new_ligands, old_keys):
            (c_id, r_name, r_id) = ok.split(".")
            lig.set_chain(c_id)
            lig.name = r_name
            lig.num_ = int(r_id)  # TODO(CJ). make this a method for the Ligand() class
            lig.residue_key = ok
            # print(lig)
            new_stru.add_chain(Chain(lig.chain(), [lig]))
        # new_stru.add(new_ligands, sort=0)
    return new_stru


def protonate_ligand(
    ligand: Ligand,
    dirname: str = ".",
    method: str = "PYBEL",
    ph: float = 7.0,
    keep_name: bool = True,
) -> Ligand:
    """Helper method that protonates a given ligand and returns the modified version."""

    _MAPPER = {
        "OPENBABEL": _protonate_ligand_OPENBABEL,
        "PYBEL": _protonate_ligand_PYBEL,
    }

    if method not in _MAPPER:
        error_msg = f"{method} is not supported for protonate_ligand()."
        error_msg += f"Allowed are '{', '.join(list(_MAPPER.keys()))}'"
        raise core.exception.UnsupportedMethod(error_msg)
    else:
        # TODO(CJ): Do we want to delete the temporary files?
        path = f"{dirname}/ligand_{ligand.name}.pdb"
        fs.safe_mkdir(dirname)
        ligand.build(path)
        base_name = fs.base_file_name(path)
        out_path = f"{dirname}/{base_name}_aH.pdb"
        return _MAPPER[method](path, ph, out_path)

def _protonate_ligand_OPENBABEL(path: str, ph: float, out_path: str) -> None:
    raise Exception(f"Method: __protonate_OPENBABEL is not implemented yet!")


def _ob_pdb_charge(pdb_path: str) -> int:
    # TODO(CJ): add tests for this function
    """
    extract net charge from openbabel exported pdb file
    """
    pdb_ls = read_pdb_lines(pdb_path)
    net_charge = 0
    for pdb_l in pdb_ls:
        if pdb_l.is_HETATM() or pdb_l.is_ATOM():
            raw: str = pdb_l.get_charge()
            raw = raw.strip()
            if not len(raw):
                continue
            charge = pdb_l.charge[::-1]
            core._LOGGER.info(
                f"Found formal charge: {pdb_l.atom_name} {charge}"
            )  # TODO make this more intuitive/make sense
            net_charge += int(charge)
    return net_charge


def _protonate_ligand_PYBEL(path: str, ph: float, out_path: str) -> Ligand:
    """Impelemntation of protonate_ligand() using the pybel method. SHOULD NOT be called directly by users."""

    def _fix_ob_output(pdb_path, out_path, ref_name_path=None) -> None:
        """
        fix atom label in pdb_pat write to out_path
        ---------
        ref_name_path: if use original atom names from pdb
        - default: None
            according to tleap output, the name could be just *counting* the element start from ' ' to number
        - : not None
            check if there're duplicated names originally, add suffix if there are.
        """
        if ref_name_path:
            ref_a_atoms = read_pdb_lines(ref_name_path)
            pdb_atoms = read_pdb_lines(pdb_path)
            ref_a_atoms = list(
                filter(lambda aa: aa.is_ATOM() or aa.is_HETATM(), ref_a_atoms)
            )
            pdb_atoms = list(
                filter(lambda aa: aa.is_ATOM() or aa.is_HETATM(), pdb_atoms)
            )
            assert len(pdb_atoms) == len(ref_a_atoms)
            for idx, (p, a) in enumerate(zip(pdb_atoms, ref_a_atoms)):
                pdb_atoms[idx].atom_name = a.atom_name
                pdb_atoms[idx].atom_id = a.atom_id
                pdb_atoms[idx].resi_id = a.resi_id
                pdb_atoms[idx].resi_name = a.resi_name
            # exit( 0 )
        fs.write_lines(out_path, list(map(lambda aa: aa.build(), pdb_atoms)))
    def fix_atom_naming( out_path ):
        lig_name = Path(out_path).stem.split('_')[1].strip()
        lines = read_pdb_lines(out_path)
        lines = list(filter(lambda ll: ll.is_ATOM() or ll.is_HETATM() or ll.line.startswith('END'), lines))                
        for ll in lines:
            if ll.line.startswith('END'):
                continue
            rawline = ll.line
            aname = rawline[76:78].strip() 
            ll.line = rawline[0:13] + f"{aname: <4}{lig_name: >3}"+rawline[20:]
        fs.write_lines(out_path, list(map(lambda ll: ll.line, lines)))
    pybel.ob.obErrorLog.SetOutputLevel(0)
    mol = next(pybel.readfile("pdb", path))
    mol.OBMol.AddHydrogens(False, True, ph)
    mol.write("pdb", out_path, overwrite=True)
    #_fix_ob_output(out_path, out_path, path)
    fix_atom_naming(out_path)
    return ligand_from_pdb(out_path, _ob_pdb_charge(out_path))

PEPTIDE_PROTONATION_METHODS = {
    "pdb2pqr" : protonate_peptide_with_pdb2pqr
}

LIGAND_PROTONATION_METHODS = {
    "pybel" : protonate_ligand_with_pybel
}

