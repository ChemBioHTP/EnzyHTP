"""Defines functions for protonating Structure(), Ligand() and MetalAtom() objects as well as raw PDB files.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-05
"""
#TODO(CJ): add more documentation
import enzy_htp.core as core
from enzy_htp.core import file_system as fs
from enzy_htp.structure import Structure, Ligand, structure_from_pdb, ligand_from_pdb, Chain

from pdb2pqr.main import main_driver as run_pdb2pqr
from pdb2pqr.main import build_main_parser as build_pdb2pqr_parser
import openbabel
import openbabel.pybel as pybel
from .pdb_line import read_pdb_lines


def check_valid_ph( ph: float ) -> None: 
    """Helper function that checks the pH is on the range: [0.00, 14.00]. Throws core.InvalidPH() if not."""
    if ph < 0 or ph > 14:
        raise core.InvalidPH(f"{ph:.2f} must be on range: [0.00,14.00]")


def protonate_pdb( pdb_path : str, pqr_path : str, ph : float = 7.0, ffout : str = "AMBER"):
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


def protonate_missing_elements(old_pdb_path : str, new_pdb_path : str, work_dir : str) -> Structure:
    """Method that compares two Structure() objects and combines missin elements from old_str to new_stru"""
    # TODO(CJ) Maybe this should be part of the structure.structure.py file
    old_stru = structure_from_pdb(old_pdb_path)
    new_stru = structure_from_pdb(new_pdb_path)
    
    if old_stru == new_stru: #CHECK(CJ): Only sequence elements.
        return new_stru
        new_stru = self.merge_structure_elements(old_stru, new_stru)

    #print(compare_structures(old_stru,new_stru))
    metal_list = old_stru.get_metals()
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
    ligand_list = old_stru.get_ligands()
    if len(ligand_list):
        print('lig-list',ligand_list)
        core._LOGGER.info(f"Merging {len(ligand_list)} ligands in old structure!")
        lig_dir = work_dir + "/ligands/"
        fs.safe_mkdir(lig_dir)
        # TODO: logging
        #print(lig_dir)
        #print("-asdgasdgasg")
        old_keys = list(map(lambda l: l.residue_key, ligand_list))
        new_ligands = list(
            map(
                lambda ll: protonate_ligand(ll, dirname=lig_dir, ph=7),
                ligand_list,
            )
        )
        print(new_ligands)
        for lig, ok in zip(new_ligands, old_keys):
            (c_id, r_name, r_id) = ok.split(".")
            lig.set_chain(c_id)
            lig.name = r_name
            lig.num_ = int(
                r_id
            )  # TODO(CJ). make this a method for the Ligand() class
            lig.residue_key = ok
            print(lig)
            new_stru.insert_chain(Chain(lig.chain(), [lig]))
        # new_stru.add(new_ligands, sort=0)
    return new_stru


def protonate_ligand(
    ligand: Ligand, dirname=".", method="PYBEL", ph=7.0, keep_name=1
) -> Ligand:
    # def protonate_ligand(cls, path, method="PYBEL", ph=7.0, keep_name=1):
    """
    Protonate the ligand from 'path' with 'method', provide out_path and net charge.
    TODO "obabel -ipdb ligand_1.pdb -opdb pdb -O ligand_1_aHt.pdb -h" can keep names, but how is it accessed by pybel
    ---------------
    method      : PYBEL (default)
                  Dimorphite (from https://durrantlab.pitt.edu/dimorphite-dl/) TODO seems better and with better python API.
                  OPENBABEL (not working if block warning output)
    ph          : 7.0 by default 
    keep_name   : if keep original atom names of ligands (default: 1)
                    - check if there're duplicated names, add suffix if are.
    """
    path = f"{dirname}/ligand_{ligand.name}.pdb"
    ligand.build(path)
    outp1_path = path[:-4] + "_badname_aH.pdb"
    out_path = path[:-4] + "_aH.pdb"
    # outm2_path = path[:-4]+'_aH.mol2'

    if method == "OPENBABEL":
        # not working if block warning output for some reason
        # openbabel.obErrorLog.SetOutputLevel(0)
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "pdb")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, path)
        mol.AddHydrogens(False, True, ph)
        obConversion.WriteFile(mol, out_path)
    if method == "PYBEL":
        pybel.ob.obErrorLog.SetOutputLevel(0)
        mol = next(pybel.readfile("pdb", path))
        mol.OBMol.AddHydrogens(False, True, ph)
        mol.write("pdb", outp1_path, overwrite=True)
        # fix atom label abd determing net charge
        if keep_name:
            _fix_ob_output(outp1_path, out_path, ref_name_path=path)
        else:
            _fix_ob_output(outp1_path, out_path)
        # determine partial charge
        # > METHOD 1<
        net_charge = _ob_pdb_charge(outp1_path)
        # > METHOD 2 <
        # mol.write('mol2', outm2_path, overwrite=True)
        # mol = next(pybel.readfile('mol2', outm2_path))
        # net_charge=0
        # for atom in mol:
        #     net_charge=net_charge+atom.formalcharge
    if method == "Dimorphite":
        pass
    return ligand_from_pdb(out_path, net_charge)


def protonation_metal_fix(self, Fix): #TODO(CJ): change to protonate_metal()
    """
    return a bool: if there's any metal center
    """
    # try once if not exist
    if self.metal_centers == []:
        self.get_metal_center()
    if self.metal_centers == []:
        print("No metal center is found. Exit Fix.")
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


def _fix_ob_output(pdb_path, out_path, ref_name_path=None):
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

#            assert len(pdb_atoms) == len(ref_a_names )
#            print('hereere')
#            pdb_ls = read_pdb_lines(pdb_path)
#            ref_resi_name = pdb_ls[0].resi_name
#            for pdb_l in pdb_ls:
#                if pdb_l.is_HETATM() or pdb_l.is_ATOM():
#                    # pybel use line order (not atom id) to assign new atom id
#                    ref_a_names.append(pdb_l.atom_name)
#        print(ref_a_names)
#        # count element in a dict
#        ele_count = {}
#        pdb_ls = read_pdb_lines(pdb_path)
#        line_count = 0
#        lines = []
#        for pdb_l in pdb_ls:
#            if not pdb_l.is_ATOM() and not pdb_l.is_HETATM():
#                continue
#
#            if ref_name_path == None:
#                ele = pdb_l.get_element()
#            else:
#                if line_count < len(ref_a_names):
#                    ele = ref_a_names[line_count]
#                else:
#                    ele = pdb_l.get_element()  # New atoms
#                pdb_l.resi_name = ref_resi_name
#                line_count += 1
#            # determine the element count
#            try:
#                # rename if more than one (add count)
#                ele_count[ele] += 1
#                pdb_l.atom_name = ele + str(ele_count[ele])
#            except KeyError:
#                ele_count[ele] = 0
#                pdb_l.atom_name = ele
#            lines.append(pdb_l.build())

# write_lines(out_path, lines)


def _ob_pdb_charge(pdb_path):
    """
    extract net charge from openbabel exported pdb file
    """

    pdb_ls = read_pdb_lines(pdb_path)
    net_charge = 0
    for pdb_l in pdb_ls:
        if pdb_l.is_HETATM() or pdb_l.is_ATOM():
            if len(pdb_l.get_charge()) != 0:
                charge = pdb_l.charge[::-1]
                core._LOGGER.info(
                    f"Found formal charge: {pdb_l.atom_name} {charge}"
                )  # TODO make this more intuitive/make sense
                net_charge += int(charge)
    return net_charge

