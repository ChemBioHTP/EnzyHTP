"""Definition for the Atom class. Meant to be the base unit of structure in package. Most of this 
functionality should NOT be directly called by users, but used within algorithm development instead.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
from __future__ import annotations
from typing import Union, Any

import numpy as np


class Atom:
    """Base unit of structure in enzy_htp. Effectively acts as dict() to store various data values 
    from each PDB lines. Designed to be used in tandem with the biopandas.pdb.PandasPdb() class.
    Each Atom() instance is created directly from each row in a PandasPDB() DataFrame(). 

    Note: The Atom class SHOULD NOT be created directly by users. It is handled by other functions in enzy_htp.
   
    Attributes:
		alt_loc : Character indicating if alternative location is present.
		atom_name : The name of the atom as a string.
		atom_number : An integer number for the atom.
		b_factor : A float representing the temperature factor (b factor).
		blank_1 : Filler space 1.
		blank_2 : Filler space 2.
		blank_3 : Filler space 3.
		blank_4 : Filler space 4.
		chain_id : Optional letter designating parent chain.
		charge : Charge of the atom.
		element_symbol : Character representing element.
		insertion : Character chode for insertion of residues.
		line_idx : Index of line from original .pdb file.
		occupancy : Float representing occupancy of the atom.
		record_name : The record type of the atom as a string (ATOM, HETATM, etc.).
		residue_name : Name of parent residue as string.
		residue_number : Number for parent residue as integer.
		segment_id : Segment identifier as a string.
		x_coord : x coordinate in cartesian 3D space.
		y_coord : y coordinate in cartesian 3D space.
		z_coord : z coordinate in cartesian 3D space.
    """

    def __init__(self, **kwargs):
        """Constructor meant to be called as Atom(**row), where row is of type pandas.Series"""
        self.alt_loc = str()
        self.atom_name = str()
        self.atom_number = int()
        self.b_factor = float()
        self.blank_1 = None
        self.blank_2 = None
        self.blank_3 = None
        self.blank_4 = None
        self.chain_id = str()
        self.charge = float()
        self.element_symbol = str()
        self.insertion = None
        self.line_idx = int()
        self.occupancy = float()
        self.record_name = str()
        self.residue_name = str()
        self.residue_number = int()
        self.segment_id = int()
        self.x_coord = float()
        self.y_coord = float()
        self.z_coord = float()

        for k, v in kwargs.items():
            setattr(self, k, v)

    def residue_key(self) -> str:
        """Creates a residue key of format "chain_id.residue_name.residue_number". Used for residue construction."""
        return f"{self.chain_id}.{self.residue_name}.{self.residue_number}"

    def distance(self, other: Atom) -> float:
        """Calculates 3D space between two atoms."""
        return np.sqrt(
            (self.x_coord - other.x_coord) ** 2
            + (self.y_coord - other.y_coord) ** 2
            + (self.z_coord - other.z_coord) ** 2
        )

    def __getitem__(self, key: str) -> Union[None, bool, int, float, str]:
        """Accessor to make Atom() callable like dctionary. If attribute is not defined, returns None """
        return self.__dict__.get(key, None)

    def __setitem__(self, key: str, value: Any) -> None:
        """Setter to make indexable like a dicT()."""
        self.__dict__[key] = value

    def to_pdb_line(
        self,
        a_id: int = None,
        r_id: int = None,
        c_id: str = None,
        ff: str = "AMBER",
        forcefield: str = "ff14SB",
    ) -> str:
        """I/O method for converting an Atom() to a PDB line. User specifies both force field type and specific force field."""
        # default
        if a_id == None:
            a_id = self.atom_number
        if r_id == None:
            r_id = self.residue_number
        if c_id == None:
            c_id = self.chain_id
        a_type = self.atom_name
        if ff == "AMBER":
            # build an amber style line here
            l_type = "{:<6}".format("ATOM")
            a_index = "{:>5d}".format(a_id)

            if forcefield == "ff14SB":
                # ff14SB by default
                # atom name TODO when deal with 2-letter element in artifical residue and ligand
                if len(self.atom_name) > 3:
                    a_name = "{:<4}".format(self.atom_name)
                else:
                    a_name = "{:<3}".format(self.atom_name)
                    a_name = " " + a_name
                r_name = "{:>3}".format(self.residue_name)
            else:
                raise Exception("Only support ff14SB atom/resiude name now")

            c_index = c_id.strip()
            if not c_index:
                c_index = " "
            r_index = "{:>4d}".format(r_id)
            r_name = r_name.strip()
            x = "{:>8.3f}".format(self.x_coord)
            y = "{:>8.3f}".format(self.y_coord)
            z = "{:>8.3f}".format(self.z_coord)
        # TODO(CJ) figure out a less janky way to do all of this
        a_type = "".join(list(filter(lambda c: not c.isdigit(), a_name))).strip()[0]
        return f"{l_type}{a_index} {a_name} {r_name} {c_index}{r_index}    {x}{y}{z}  1.00  0.00           {a_type}  "

    #

    def build_oniom(self, layer, chrg=None, cnt_info: list = None, if_lig=0, if_sol=0):
        """
        build line for oniom. Use element name for ligand atoms since they are mostly in QM regions.
        Gaussian use *ff96* which the atom type is corresponding to all_amino94.lib and ion94.lib in Amber distribution
        ---------
        chrg    : charge of current atom
        layer   : layer where atom in
        cnt_info: (low atom only) when current atom was a boundary atom. Provide [cnt_atom_ele, cnt_atom_label, cnt_atom_id]
        """
        cnt_flag = ""
        self.get_ele()
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

        # label: deal with N/C terminal and ligand
        if if_lig:
            G16_label = self.ele
        if if_sol:
            G16_label = G16_label_map[self.parent.name][self.name]
        if not if_sol and not if_lig:
            r = self.resi
            r1 = self.resi.chain.residues[0]
            rm1 = self.resi.chain.residues[-1]
            if r == r1 and self.name in ["H1", "H2", "H3"]:
                G16_label = "H"
            else:
                if r == rm1 and self.name == "OXT":
                    G16_label = "O2"
                else:
                    G16_label = G16_label_map[self.parent.name][self.name]

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

        return f"{atom_label} {fz_flag}   {x} {y} {z} {ly_flag}{cnt_flag}{line_feed}"
