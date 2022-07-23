"""Defines functions for protonating Structure(), Ligand() and MetalAtom() objects as well as raw PDB files.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-05
"""
# TODO(CJ): add more documentation
from pathlib import Path

import enzy_htp.core as core
from enzy_htp.core import file_system as fs
from enzy_htp.structure import (
    Structure,
    Ligand,
    structure_from_pdb,
    ligand_from_pdb,
    Chain,
)

from pdb2pqr.main import main_driver as run_pdb2pqr
from pdb2pqr.main import build_main_parser as build_pdb2pqr_parser
import openbabel
import openbabel.pybel as pybel
from .pdb_line import read_pdb_lines

# TODO(CJ): this probably needs to go to core so that structure.Structure can use it
def check_valid_ph(ph: float) -> None:
    """Helper function that checks the pH is on the range: [0.00, 14.00]. Throws core.InvalidPH() if not."""
    if ph < 0 or ph > 14:
        raise core.InvalidPH(f"{ph:.2f} must be on range: [0.00,14.00]")


def protonate_pdb(
    pdb_path: str, pqr_path: str, ph: float = 7.0, ffout: str = "AMBER"
) -> None:
    """Runs PDB2PQR on a specified pdb file and saves it to the specified pqr path. This preparation step
    makes use of [PDB2PQR](https://www.poissonboltzmann.org/) via the pdb2pqr python [package](https://pdb2pqr.readthedocs.io/en/latest/).
    Adds in missing atoms and finds the protonation state of the pdb file.
    """
    # TODO(CJ): check if ffout is valid.
    # TODO(CJ): maybe improve the documentation here?
    check_valid_ph(ph)
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

    # print(compare_structures(old_stru,new_stru))
    metal_list = old_stru.metals
    if len(metal_list):
        core._LOGGER.info(f"Merging {len(metal_list)} metal centers in old structure!")
        core._LOGGER.info(f"Adding metal centers to new structure...")
        new_stru.add(metal_list, sort=0)
        core._LOGGER.info(f"Metal centers added!")
        # fix metal environment
        core._LOGGER.info(f"Protonating newly added metals...")
        new_stru.protonation_metal_fxi(Fix=1)
        core._LOGGER.info(f"Protonation complete!")

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
            new_stru.insert_chain(Chain(lig.chain(), [lig]))
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


def protonation_metal_fix(self, Fix):  # TODO(CJ): change to protonate_metal()
    """
    return a bool: if there's any metal center
    """
    # try once if not exist
    if self.metal_centers == []:
        self.get_metal_center()
    if self.metal_centers == []:
        print("No meNontal center is found. Exit Fix.")
        return False

    # start fix
    # get donor atoms and residues
    for metal in self.metal_centers:
        metal.get_donor_residue(method="INC")

        if Fix == 1:
            metal._metal_fix_1()

        if Fix == 2:
            metal._metal_fix_2()

        if Fix == 3:
            metal._metal_fix_3()
    return True


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
