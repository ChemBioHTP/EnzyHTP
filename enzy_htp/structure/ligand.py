from .atom import Atom
from typing import List
from .residue import Residue
from .pdb_line import PDBLine, read_pdb_lines
from ..core import get_file_ext, write_lines, UnsupportedFileType, lines_from_file, _LOGGER

import openbabel
import openbabel.pybel as pybel


class Ligand(Residue):
    """
    -------------
    initilize from
    PDB:        Residue.fromPDB(atom_input, input_type='PDB_line' or 'line_str' or 'file' or 'path')
    raw data:   Residue(atom_name, coord)
    Residue:    Ligand.fromResidue(atom_obj)
    -------------
    id
    name
    atoms = [atom, ...]
    parent # the whole stru
    -------------
    Method
    -------------
    set_parent
    -------------
    """

    def __init__(self, residue_key: str, atoms: List[Atom], **kwargs):
        self.net_charge = kwargs.get('net_charge', None)
        Residue.__init__(self, 'A.B.10', atoms)
#
#    @classmethod
#
#    @classmethod
#    def fromPDB(
#        cls,
#        resi_input,
#        resi_id=None,
#        resi_name=None,
#        net_charge=None,
#        input_type="PDB_line",
#    ):
#        """
#        generate resi from PDB. Require 'ATOM' and 'HETATM' lines.
#        ---------
#        resi_input = PDB_line (or line_str or file or path)
#        resi_id : int (use the number in the line by default // support customize)
#        net_charge : user assigned net charge for further use
#        Use PDB_line in the list to init each atom
#        """
#
#        # adapt general input // converge to a list of PDB_line (resi_lines)
#        if input_type == "path":
#            f = open(resi_input)
#            resi_lines = PDB_line.fromlines(f.read())
#            f.close()
#        if input_type == "file":
#            resi_lines = PDB_line.fromlines(resi_input.read())
#        if input_type == "line_str":
#            resi_lines = PDB_line.fromlines(resi_input)
#        if input_type == "PDB_line":
#            resi_lines = resi_input
#
#        # clean lines
#        for i in range(len(resi_lines) - 1, -1, -1):
#            if (
#                resi_lines[i].line_type != "ATOM"
#                and resi_lines[i].line_type != "HETATM"
#            ):
#                if Config.debug > 1:
#                    print("Residue.fromPDB: delete error line in input.")
#                del resi_lines[i]
#
#        # Default resi_id
#        if resi_id is None:
#            resi_id = resi_lines[0].resi_id
#        # get name from first line
#        if resi_name is None:
#            resi_name = resi_lines[0].resi_name
#        # get child atoms
#        atoms = []
#        for pdb_l in resi_lines:
#            atoms.append(Atom.fromPDB(pdb_l))
#
#        return cls(atoms, resi_id, resi_name, net_charge=net_charge)
#
#    def sort(self):
#        pass
#

    def get_net_charge(self) -> int:
        return self.net_charge

    def build(self, out_path: str) -> None:
        """
        build ligand file(out_path) with ft format
        """
        ext = get_file_ext(out_path).lower()
        if ext != '.pdb':
            raise UnsupportedFileType(out_path)
        lines = list(
            map(lambda pr: pr[1].build(a_id=pr[0] + 1, c_id=" "),
                enumerate(self.atoms))) + ["TER", "END"]
        write_lines(out_path, lines)
        return
        with open(out_path, "w") as of:
            a_id = 0
            for atom in self:
                a_id = a_id + 1
                line = atom.build(a_id=a_id, c_id=" ")
                of.write(line)
            of.write("TER" + line_feed + "END" + line_feed)


#    def get_net_charge(self, method="PYBEL", ph=7.0, o_dir="."):
#        """
#        get net charge for the ligand
#        -------
#        method   : PYBEL (default) use UNITY_ATOM_ATTR info from openbabel mol2
#        """
#        # build file
#        mkdir(o_dir + "/cache")
#        temp_path = o_dir + "/cache/ligand_temp.pdb"
#        temp_pdb2_path = o_dir + "/cache/ligand_temp2.pdb"
#        temp_pdb3_path = o_dir + "/cache/ligand_temp3.pdb"
#        temp_m2_path = o_dir + "/cache/ligand_temp.mol2"
#        self.build(temp_path)
#        # charge
#        if method == "PYBEL":
#            pybel.ob.obErrorLog.SetOutputLevel(0)
#            # remove H (or the )
#            mol = next(pybel.readfile("pdb", temp_path))
#            mol.removeh()
#            mol.write("pdb", temp_pdb2_path, overwrite=True)
#            # clean connectivity
#            os.system("cat " + temp_pdb2_path + " |grep 'ATOM' > " + temp_pdb3_path)
#            # add H and result net charge
#            mol = next(pybel.readfile("pdb", temp_pdb3_path))
#            mol.OBMol.AddHydrogens(False, True, ph)
#            mol.write("mol2", temp_m2_path, overwrite=True)
#            mol = next(pybel.readfile("mol2", temp_m2_path))
#            net_charge = 0
#            for atom in mol:
#                net_charge = net_charge + atom.formalcharge
#        self.net_charge = net_charge
#        return net_charge
#


def residue_to_ligand(ptr: Residue) -> Ligand:
    """
    generate from a Residue object. copy data
    """
    return Ligand(ptr.name, ptr.atoms)


def protonate_ligand(ligand: Ligand,
                     dirname='.',
                     method="PYBEL",
                     ph=7.0,
                     keep_name=1) -> Ligand:
    #def protonate_ligand(cls, path, method="PYBEL", ph=7.0, keep_name=1):
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
    return out_path, net_charge


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
        ref_a_names = []
        pdb_ls = read_pdb_lines(pdb_path)
        ref_resi_name = pdb_ls[0].resi_name
        for pdb_l in pdb_ls:
            if pdb_l.line_type == "HETATM" or pdb_l.line_type == "ATOM":
                # pybel use line order (not atom id) to assign new atom id
                ref_a_names.append(pdb_l.atom_name)

    # count element in a dict
    ele_count = {}
    pdb_ls = read_pdb_lines(pdb_path)
    line_count = 0
    lines = []
    for pdb_l in pdb_ls:
        if not pdb_l.is_ATOM() and not pdb_l.is_HETATM():
            continue

        if ref_name_path == None:
            ele = pdb_l.get_element()
        else:
            if line_count < len(ref_a_names):
                ele = ref_a_names[line_count]
            else:
                ele = pdb_l.get_element()  # New atoms
            pdb_l.resi_name = ref_resi_name
            line_count += 1
        # determine the element count
        try:
            # rename if more than one (add count)
            ele_count[ele] += 1
            pdb_l.atom_name = ele + str(ele_count[ele])
        except KeyError:
            ele_count[ele] = 0
            pdb_l.atom_name = ele
        lines.append(pdb_l.build())

    write_lines(out_path, lines)


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
                _LOGGER.info(f"Found formal charge: {pdb_l.atom_name} {charge}"
                            )  # TODO make this more intuitive/make sense
                net_charge += int(charge)
    return net_charge
