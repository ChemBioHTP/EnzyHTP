"""Definition for the MetalUnit class. This is a specialization of the Residue() obejct that is
stored alongside other Residue() children types in the Chain() class.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""
from __future__ import annotations

from copy import deepcopy
import sys
import numpy as np
from typing import List, Union
import enzy_htp.chemical as chem
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Residue, Atom

class MetalUnit(Residue):
    """Class representing a residue unit that contains a meetal atom in a protein/enzyme
    structure or system. Only a single atom is allowed since atoms are assumed to have connectivity
    in the residue

    Attributes:
        (same as residue)
    Derived properties:
        (same as residue)
        atom : the only one atom in this MetalUnit unit
        atom_name : the name of the atom of the MetalUnit unit
    """

    def __init__(self, residue_idx: int, residue_name: str, atoms: List[Atom], parent=None):
        """Constructor for MetalUnit. Identical to Residue() constructor."""
        Residue.__init__(self, residue_idx, residue_name, atoms, parent)
        self.set_rtype(chem.ResidueType.METAL)

    # === Getter-Attr (ref) ===
    @property
    def atom(self):
        """Getter for the only one atom in this MetalUnit unit"""
        return self._atoms[0]
    @atom.setter
    def atom(self, val):
        self._atoms[0] = val

    @property
    def atom_name(self):
        """Getter for the only one atom in this MetalUnit unit"""
        return self._atom[0].name
    @atom_name.setter
    def atom_name(self, val):
        self._atom[0].name = val

    # TODO below to change

    def location(self) -> np.array: #@shaoqz: maybe unify as coord etc. TODO
        """Gets the location of the MetalUnit."""
        return self._atoms[0]

    # === Getter-Prop (cpy/new) ===
    def get_radii(self, method: str = "ionic") -> Union[float, None]:
        """Gets the atomic radii for the metal using specified method.
        Allowed values are 'ionic' or 'vdw' for ionic and van-der waals, respectively.
        Returns a float or None if the metal does not have a distance.
        """
        return chem.get_metal_radii(self.atom_name_, method)

    def clone(self) -> MetalUnit:
        """Creates deepcopy of self."""
        return deepcopy(self)

    # === Checker === 
    def is_metal(self) -> bool:
        """Checks if Residue is a metal. Always returns True for this specialization."""
        return True

    def is_metal_center(self) -> bool:
        """Checks if Residue is a metal center. Always returns True for this specialization."""
        return True #@shaoqz: @imp2 not all metal is metalcenter. e.g.: Na+ and K+ is not considered a metal coordinating center.

    def is_canonical(self) -> bool:
        """Checks if the Residue is a canonical. Always returns False for this specialization."""
        return False

    # === Editor ===
    # === Special ===
    def __getitem__(self, key: int):
        """
        Metalatom_obj[0]: self
        -----
        return self when iter. Allow metalatom to be iterate like a residue that return a atom level obj.
        """
        if key == 0:
            return self
        if key > 0:
            raise StopIteration

    #region === TODO/TOMOVE ===
    # fix related
    def get_donor_atom(self, method: str = "INC", check_radius=4.0): # @nu
        # TODO(CJ): move this higher up in the structural hierarchy.
        """
        Get coordinated donor atom for a metal center.
        1. check all atoms by type, consider those in the "donor_map"
        (only those atoms from residues are considered. So atoms from ligand and ion will be ignored)
        2. check if atoms are within the check_radius
        3. check distance for every atoms left.
        -----------------------
        Base on a certain type of radius of this metal:
        method = INC ------ {ionic radius} for both metal and donor atom
               = VDW ------ {Van Der Waals radius} for both metal and donor atom
        check_radius ------ the radius that contains all considered potential donor atoms. set for reducing the complexity.
        -----------------------
        save found atom list to self.donor_atoms
        """

        # find radius for matal
        if method == "INC":
            R_m = Ionic_radius_map[self.ele]
        if method == "VDW":
            R_m = VDW_radius_map[self.ele]

        # get target with in check_radius (default: 4A)
        coord_m = np.array(self.coord)
        protein_atoms = self.parent.get_all_protein_atom()
        for atom in protein_atoms:

            # only check donor atom (by atom name)
            if atom.name in Donor_atom_list[atom.ff]:

                # cut by check_radius
                dist = np.linalg.norm(np.array(atom.coord) - coord_m)
                if dist <= check_radius:
                    # determine coordination
                    atom.get_ele()
                    if method == "INC":
                        R_d = Ionic_radius_map[atom.ele]
                    if method == "VDW":
                        R_d = VDW_radius_map[atom.ele]

                    if dist <= (R_d + R_m):
                        self.donor_atoms.append(atom)
                        if Config.debug > 1:
                            print(
                                "Metalatom.get_donor_atom: "
                                + self.name
                                + " find donor atom:"
                                + atom.resi.name
                                + " "
                                + str(atom.resi.id)
                                + " "
                                + atom.name
                            )

    def get_donor_residue(self, method="INC"): # @nu
        """
        get donor residue based on donor atoms
        """
        self.donor_resi = []
        self.get_donor_atom(method=method)
        # add d_atom and save donor_resi
        for atom in self.donor_atoms:
            resi = atom.resi
            resi.d_atom = atom
            resi.a_metal = self
            self.donor_resi.append(resi)

        # warn if more than one atom are from a same residue
        for index in range(len(self.donor_resi)):
            for index2 in range(len(self.donor_resi)):
                if index2 > index:
                    if self.donor_resi[index2].id == self.donor_resi[index].id:
                        print(
                            "\033[1;31;40m!WARNING! found more than 1 donor atom from residue: "
                            + self.donor_resi[index].name
                            + str(self.donor_resi[index].id)
                            + "\033[m"
                        )

    def _metal_fix_1(self): # @nu
        """
        Fix1: deprotonate all donor (rotate those with tight H, like Ser)
            - warn if uncommon donor residue shows up (like Ser)
        """
        for resi in self.donor_resi:
            if resi.ifDeProton():
                resi.deprotonate(resi.d_atom)
            else:
                if resi.name not in NoProton_list:
                    print(
                        "!WARNING!: uncommon donor residue -- "
                        + resi.chain.id
                        + " "
                        + resi.name
                        + str(resi.id)
                    )
                    # resi.rot_proton(resi.d_atom)

    def _metal_fix_2(self): # @nu
        """
        Fix2: rotate if there're still lone pair left
        """
        for resi in self.donor_resi:
            resi.rot_proton(resi.d_atom)

    def _metal_fix_3(self): # @nu
        """
        Fix3: run pka calculate containing ion (maybe pypka) and run Fix2 based on the result
        waiting for response
        """
        pass

    def build(self, a_id=None, r_id=None, c_id=None, ff="AMBER", forcefield="ff14SB"): # @nu
        """
        generate an metal atom output line. End with LF
        return a line str
        --------
        use self.id if not assigned

        """
        # default
        if a_id == None:
            a_id = self.id
        if r_id == None:
            r_id = self.resi.id
        if c_id == None:
            print("WARNING: the metal atom may need a chain id in build()!")
            c_id = " "

        if ff == "AMBER":
            # build an amber style line here
            l_type = "{:<6}".format("ATOM")
            a_index = "{:>5d}".format(a_id)

            if forcefield == "ff14SB":
                # ff14SB by default
                a_name = "{:>2}".format(self.name)
                r_name = "{:<3}".format(self.resi_name)
            else:
                raise Exception("Only support ff14SB atom/resiude name now")

            c_index = c_id
            r_index = "{:>4d}".format(r_id)
            x = "{:>8.3f}".format(self.coord[0])
            y = "{:>8.3f}".format(self.coord[1])
            z = "{:>8.3f}".format(self.coord[2])

        # example: ATOM   5350  HB2 PRO   347      32.611  15.301  24.034  1.00  0.00
        line = (
            l_type
            + a_index
            + " "
            + a_name
            + "   "
            + r_name
            + " "
            + c_index
            + r_index
            + "    "
            + x
            + y
            + z
            + "  1.00  0.00"
            + line_feed
        )

        return line

    def build_oniom(self, layer, chrg=None, cnt_info: list = None): # @nu
        """
        build a metal line for oniom.
        Gaussian use *ff96* which the atom type is corresponding to all_amino94.lib and ion94.lib in Amber distribution
        For metals that do not exist in ion94.lib. Custom a type and get parms from TIP3P lib.
        ---------
        chrg    : charge of current metal atom
        layer   : layer where atom in
        cnt_info: (low atom only) when current atom was a boundary atom. Provide [cnt_atom_ele, cnt_atom_label, cnt_atom_id]
        """
        cnt_flag = ""
        if layer not in ["h", "l"]:
            raise Exception('build_oniom: please use: "h" or "l" for layer')
        if layer == "h":
            fz_flag = "0"
            ly_flag = "H"
        if layer == "l":
            fz_flag = "-1"
            ly_flag = "L"
            if cnt_info != None:
                cnt_flag = (
                    " " + cnt_info[0] + "-" + cnt_info[1] + " " + str(cnt_info[2])
                )

        # label
        if self.resi_name in G16_label_map.keys():
            G16_label = G16_label_map[self.resi_name][self.name]
        else:
            if Config.debug >= 1:
                print("Metal: " + self.name + " not in build-in atom type of ff96.")
                print(
                    "Use parameters and atom types from TIP3P (frcmod.ionsjc_tip3p & frcmod.ions234lm_126_tip3p)"
                )
            G16_label = self.resi_name.strip("+-") + "0"
            self.parm = tip3p_metal_map[self.resi_name][self.name]
            self.parm[0] = G16_label

        # chrg
        if chrg == None:
            try:
                chrg = self.charge  # from prmtop
            except NameError:
                raise Exception(
                    "You need to at least provide a charge or use get_atom_charge to get one from prmtop file."
                )

        atom_label = "{:<16}".format(
            " " + self.ele + "-" + G16_label + "-" + str(round(chrg, 6))
        )
        fz_flag = "{:>2}".format(fz_flag)
        x = "{:<14.8f}".format(self.coord[0])
        y = "{:<14.8f}".format(self.coord[1])
        z = "{:<14.8f}".format(self.coord[2])

        line = (
            atom_label
            + " "
            + fz_flag
            + "   "
            + x
            + " "
            + y
            + " "
            + z
            + " "
            + ly_flag
            + cnt_flag
            + line_feed
        )

        return line

    def get_pdb_lines(self, _) -> List[str]:  #@shaoqz: go to IO
        """Method that gets the PDB lines for all of the Atom() objects."""
        result: List[str] = Residue.get_pdb_lines(self, False)
        for idx, line in enumerate(result):
            result[idx] = "HETATM" + line[6:]
        return result

    #endregion

def residue_to_metal(residue: Residue) -> MetalUnit:
    '''Convenience function that converts Residue() to MetalUnit() object.'''
    if len(residue.atoms) > 1:
        _LOGGER.error(f'Found more than 1 atom in a metal residue unit: {residue.idx} {residue.name}')
        sys.exit(1)
    return MetalUnit(residue.idx, residue.name, residue.atoms, residue.parent)
