"""Definition for the MetalUnit class. This is a specialization of the Residue() obejct that is
stored alongside other Residue() children types in the Chain() class.
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""
from __future__ import annotations

from copy import deepcopy
import itertools
import sys
from typing import List, Tuple, Union, Dict

from .noncanonical_base import NonCanonicalBase
import enzy_htp.chemical as chem
from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Residue, Atom


class MetalUnit(NonCanonicalBase):
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

    def __init__(self, residue_idx: int, residue_name: str, atoms: List[Atom], parent=None, **kwargs):
        """Constructor for MetalUnit. Identical to Residue() constructor."""
        NonCanonicalBase.__init__(self, residue_idx, residue_name, atoms, parent, **kwargs)
        self.rtype = chem.ResidueType.METAL

    #region === Getter-Attr ===
    @property
    def atom(self) -> Atom:
        """Getter for the only one atom in this MetalUnit unit"""
        return self._atoms[0]

    @atom.setter
    def atom(self, val: Atom):
        self._atoms[0] = val

    @property
    def atom_name(self) -> str:
        """Getter for the only one atom in this MetalUnit unit"""
        return self._atoms[0].name

    @atom_name.setter
    def atom_name(self, val: str):
        self._atoms[0].name = val

    #endregion

    #region === Getter-Prop ===
    @property
    def element(self) -> str:
        """determine the element name of current metalunit"""
        return chem.metal.METAL_MAPPER[self.name]

    @property
    def coord(self) -> Tuple[int, int, int]:
        """return the coordinate of the atom in the unit"""
        return self.atom.coord

    def get_donor_mapper(self, method: str = "ionic", check_radius: int = 4.0) -> Dict[Residue, List[Atom]]:
        """
        get a mapper of {donor_residue : donor_atom_in_this_residue, ...} for a metal center
        """
        donor_mapper = {}
        donor_atoms = self.get_donor_atoms(method, check_radius)
        atoms_groupby_residue = itertools.groupby(
            donor_atoms, lambda a: a.residue)  # be careful if you want to turn this into a list. groupitem may lose
        for i, j in atoms_groupby_residue:
            donor_mapper[i] = list(j)
        # san check
        for k, v in donor_mapper.items():
            if len(v) > 1:
                _LOGGER.warning(f"More than 1 donor atom in residue {k} to center {self}: {v}")
        return donor_mapper

    def get_donor_atoms(self, method: str = "ionic", check_radius=4.0) -> List[Atom]:
        """
        Get coordinated donor atom for a metal center.
        1. check all atoms by type, consider those in the "donor_map"
        (only those atoms from residues are considered. So atoms from ligand and ion will be ignored)
        2. check if atoms are within the check_radius
        3. check distance for every atoms left.
        Args:
            method: the method to determine radius of the metal and the donor atom.
            (current available keywords)
                ionic: (ionic radius) for both metal and donor atom
                vdw: (Van Der Waals radius) for both metal and donor atom
            check_radius: the radius that contains all considered potential donor atoms. set for reducing the complexity.
        Returns:
            a list of found donor atoms (in their reference)
        """
        result = []
        # find radius for matal
        r_metal = self.atom.radius(method)
        # get target with in check_radius (default: 4A)
        range_atoms = self.chain.parent.find_atoms_in_range(self.coord, range_distance=check_radius)
        atom: Atom
        for atom in range_atoms:
            # only check donor atom (by atom name)
            if atom.is_donor_atom():
                r_donor = atom.radius(method)
                if atom.distance_to(self.atom) - (r_metal + r_donor) <= 0.01:  # TODO refine this approximation
                    result.append(atom)
                    _LOGGER.info(f"found donor atom of {self}: {atom}")
        return result

    def get_donor_residue(self, method: str = "ionic", check_radius: int = 4.0):
        """
        get donor residue based on donor atoms. See get_donor_mapper for detail
        """
        return list(self.get_donor_mapper(method, check_radius).keys())

    def clone(self) -> MetalUnit:
        """Creates deepcopy of self."""
        return deepcopy(self)

    def radius(self, method: str = "ionic") -> float:
        return self.atom.radius(method)
    #endregion

    #region === Checker ===
    def is_metal_center(self):
        """determine if current metal is a coordination center.
        (place holder for override)
        TODO dispatch and make more consistant way to determine for 'boundary' metals like Mg2+"""
        return self.name in chem.METAL_CENTER_MAP

    def is_canonical(self) -> bool:
        """Checks if the Residue is a canonical. Always returns False for this specialization."""
        return False

    def is_metal(self) -> bool:
        """Checks if Residue is a metal. Always returns True for this specialization."""
        return True
    #endregion

    #region === Special ===
    def __str__(self) -> str:
        return f"MetalUnit({self._idx}, {self._name}, atom:{len(self._atoms)}, {self._parent})"

    #endregion

    #region === TODO/TOMOVE ===
    # fix related

    def _metal_fix_1(self):  # @nu
        """
        Fix1: deprotonate all donor (rotate those with tight H, like Ser)
            - warn if uncommon donor residue shows up (like Ser)
        """
        for resi in self.donor_resi:
            if resi.ifDeProton():
                resi.deprotonate(resi.d_atom)
            else:
                if resi.name not in NoProton_list:
                    print("!WARNING!: uncommon donor residue -- " + resi.chain.id + " " + resi.name + str(resi.id))
                    # resi.rot_proton(resi.d_atom)

    def _metal_fix_2(self):  # @nu
        """
        Fix2: rotate if there"re still lone pair left
        """
        for resi in self.donor_resi:
            resi.rot_proton(resi.d_atom)

    def _metal_fix_3(self):  # @nu
        """
        Fix3: run pka calculate containing ion (maybe pypka) and run Fix2 based on the result
        waiting for response
        """
        pass

    # def build_oniom(self, layer, chrg=None, cnt_info: list = None):  # @nu
    #     """
    #     build a metal line for oniom.
    #     Gaussian use *ff96* which the atom type is corresponding to all_amino94.lib and ion94.lib in Amber distribution
    #     For metals that do not exist in ion94.lib. Custom a type and get parms from TIP3P lib.
    #     ---------
    #     chrg    : charge of current metal atom
    #     layer   : layer where atom in
    #     cnt_info: (low atom only) when current atom was a boundary atom. Provide [cnt_atom_ele, cnt_atom_label, cnt_atom_id]
    #     """
    #     cnt_flag = ""
    #     if layer not in ["h", "l"]:
    #         raise Exception("build_oniom: please use: 'h' or 'l' for layer")
    #     if layer == "h":
    #         fz_flag = "0"
    #         ly_flag = "H"
    #     if layer == "l":
    #         fz_flag = "-1"
    #         ly_flag = "L"
    #         if cnt_info != None:
    #             cnt_flag = (" " + cnt_info[0] + "-" + cnt_info[1] + " " + str(cnt_info[2]))

    #     # label
    #     if self.resi_name in G16_label_map.keys():
    #         G16_label = G16_label_map[self.resi_name][self.name]
    #     else:
    #         if Config.debug >= 1:
    #             print("Metal: " + self.name + " not in build-in atom type of ff96.")
    #             print("Use parameters and atom types from TIP3P (frcmod.ionsjc_tip3p & frcmod.ions234lm_126_tip3p)")
    #         G16_label = self.resi_name.strip("+-") + "0"
    #         self.parm = tip3p_metal_map[self.resi_name][self.name]
    #         self.parm[0] = G16_label

    #     # chrg
    #     if chrg == None:
    #         try:
    #             chrg = self.charge  # from prmtop
    #         except NameError:
    #             raise Exception("You need to at least provide a charge or use get_atom_charge to get one from prmtop file.")

    #     atom_label = "{:<16}".format(" " + self.ele + "-" + G16_label + "-" + str(round(chrg, 6)))
    #     fz_flag = "{:>2}".format(fz_flag)
    #     x = "{:<14.8f}".format(self.coord[0])
    #     y = "{:<14.8f}".format(self.coord[1])
    #     z = "{:<14.8f}".format(self.coord[2])

    #     line = (atom_label + " " + fz_flag + "   " + x + " " + y + " " + z + " " + ly_flag + cnt_flag + line_feed)

    #     return line

    #endregion
    @property
    def bonds(self):
        return list()


   
def get_metal(metal_name:str, charge:int=None) -> MetalUnit:
    #@TODO(CJ): add some checks in here
    return MetalUnit(
            1,
            metal_name,
            atoms=[Atom(
                name=metal_name,
                coord=(0.0, 0.0, 0.0),
                charge=charge,
                idx=1
                )]
            )
                    
def residue_to_metal(residue: Residue) -> MetalUnit:
    """Convenience function that converts Residue() to MetalUnit() object."""
    if len(residue.atoms) > 1:
        _LOGGER.error(f"Found more than 1 atom in a metal residue unit: {residue.idx} {residue.name}")
        sys.exit(1)
    return MetalUnit(residue.idx, residue.name, residue.atoms, residue.parent)
