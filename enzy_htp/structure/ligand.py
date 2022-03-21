import numpy as np

from .atom import Atom
from typing import List
from .residue import Residue
#from enzy_htp.preparation import PDBLine, read_pdb_lines
from ..core import (
    get_file_ext,
    write_lines,
    UnsupportedFileType,
    lines_from_file,
    _LOGGER,
)

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
        self.net_charge = kwargs.get("net_charge", None)
        Residue.__init__(self, "A.B.10", atoms)
    #
    #    @classmethod
    #
    #    @classmethod
    #
    #    def sort(self):
    #        pass
    #

    def set_residue_number( self, num : int ) -> None:
        for idx, aa in enumerate( self.atoms ):
            self.atoms[idx].residue_number = num

    def get_net_charge(self) -> int:
        return self.net_charge

    def build(self, out_path: str) -> None:
        """
        build ligand file(out_path) with ft format
        """
        ext = get_file_ext(out_path).lower()
        if ext != ".pdb":
            raise UnsupportedFileType(out_path)
        lines = list(
            map(lambda pr: pr[1].to_pdb_line(a_id=pr[0] + 1, c_id=" "), enumerate(self.atoms))
        ) + ["TER", "END"]
        write_lines(out_path, lines)


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

