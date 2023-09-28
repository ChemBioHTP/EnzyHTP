"""Definition for the Atom class. Meant to be the base unit of structural information: coordinate and
topology (atom type & connectivity). Serve solely for storing and accessing data.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
from __future__ import annotations
import re
import sys
from typing import Any, List, Tuple
from plum import dispatch

import numpy as np
import pandas as pd

import enzy_htp.chemical as chem
from enzy_htp.core.exception import ResidueDontHaveAtom
import enzy_htp.core.math_helper as mh
from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode
from enzy_htp.core.logger import _LOGGER


class Atom(DoubleLinkedNode):
    """Base unit of structure in enzy_htp. Store coordinate and topology information
    in a Structure() object. Serve solely for storing and accessing data.
    * not specific to any file format that store a structure

    Attributes:
        (nessessary)
        name : The name of the atom as a string. often refer to a specific connectivity.
                TODO the name should be decouple with specific parsing logic. Current method is
                     use names in PDB format and covert every other format into this.
        coord: (x,y,z) for cartesian coordinate of the atom.
        parent/residue: the parent residue that this atom belongs to.
        (optional)
        idx: An integer number for the atomic index.
        b_factor : A float representing the temperature factor (b factor).
        charge : Charge of the atom.
        element : Character representing element.
        connect : the connectivity of this atom. (a list of reference of connected Atom objs)
    """

    def __init__(self, ds: pd.Series, parent=None):
        """Constructor of Atom(), where ds is of type pandas.Series"""
        # nessessary
        self._name = ds["atom_name"].strip()
        self._coord = (ds["x_coord"], ds["y_coord"], ds["z_coord"])
        self.set_parent(parent)
        self.set_ghost_children()
        # optional
        self._idx = None
        self._b_factor = None
        self._element = None
        self._charge = None
        self._connect = None
        ds_keys = ds.keys()
        if "atom_number" in ds_keys and not np.isnan(ds["atom_number"]):
            self._idx = ds["atom_number"]
        if "b_factor" in ds_keys and not np.isnan(ds["b_factor"]):
            self._b_factor = ds["b_factor"]
        if "element_symbol" in ds_keys and ds["element_symbol"].strip() != "":
            self._element = ds["element_symbol"].strip()
        if "charge" in ds_keys and not np.isnan(ds["charge"]):
            self._charge = float(ds["charge"])

    #region === Getter-Attr (ref) ===
    @property
    def name(self) -> str:
        """Getter for the Atom()'s name."""
        return self._name

    @name.setter
    def name(self, val):
        self._name = val

    @property
    def idx(self) -> int:
        """Getter for the Atom()'s index."""
        if self._idx is None:
            _LOGGER.warning(
                f"{sys._getframe().f_back.f_back.f_code.co_name}:{sys._getframe().f_back.f_code.co_name} is using index of an non-indexing atom: {self}"
            )  # pylint: disable=logging-fstring-interpolation,line-too-long
        return self._idx

    @idx.setter
    def idx(self, val):
        self._idx = val

    @property
    def coord(self) -> Tuple[float, float, float]:
        """Getter for coordinate of current atom"""
        return self._coord

    @coord.setter
    def coord(self, val):
        self._coord = val

    @property
    def residue(self):
        """synonym for parent"""
        return self.get_parent()

    @residue.setter
    def residue(self, val):
        self.set_parent(val)

    @property
    def b_factor(self):
        """getter for b_factor"""
        return self._b_factor

    @b_factor.setter
    def b_factor(self, val):
        self._b_factor = val

    @property
    def element(self):
        """getter for _element"""
        if self._element is None:
            if self.name in chem.residue.RESIDUE_ELEMENT_MAP["Amber"].keys():
                return chem.residue.RESIDUE_ELEMENT_MAP["Amber"][self.name]
            elif self.parent.is_metal():
                return self.parent.element
            else:
                # case: in ligand atoms are named like this H1
                return re.match("^[A-Z][a-z]?", self.name).group()
        elif self.parent.is_metal():  # in pdb some metal's element name is wrong
            return self.parent.element
        return self._element

    @element.setter
    def element(self, val):
        self._element = val

    @property
    def charge(self):
        """getter for _charge"""
        return self._charge

    @charge.setter
    def charge(self, val):
        self._charge = val

    @property
    def connect(self) -> List[Atom, str]:
        """getter for _connect, the list for (atoms, bond_type) it connects"""
        if not self.is_connected():
            _LOGGER.warning(f"There are no connection info for {self}, Initiating it.")
            self.init_connect_in_caa()  # TODO(qz): make this also works for non-caa
        return self._connect

    @connect.setter
    def connect(self, val: List[Atom, str]):
        """setter of connect is nessessary for residue level objects to set connectivity for atoms
        raise an error if called by non allowed objects.
        raise an error if the assigned val doesn't align with the format"""

        # 1. check for caller
        allowed_caller_class = ["Residue", "Atom"] # add allowed class here, prob also the adaptor class for RDKit
        # determine caller class
        caller_locals = sys._getframe(1).f_locals
        if "self" in caller_locals:
            caller_class = caller_locals["self"].__class__.__name__
            if caller_class not in allowed_caller_class:
                _LOGGER.error(f"only calling from methods from {allowed_caller_class} is allowed")
                sys.exit(1)
        
        # 2. check for data format
        type(self)._check_assigning_connect_data_type(val)

        # 3. assign the value if all passed
        self._connect = val

    @staticmethod
    def _check_assigning_connect_data_type(val: Any) -> bool:
        """Function used by the setter of connect. Check if the data format of val meetings the definiation.
        The format is defined in this function.
        Raise an error if the format is not met.
        Format: [[Atom(), "bond_order_info"], ...]"""
        check_pass = 1
        if isinstance(val, list):
            for i in val:
                if isinstance(i, list):
                    if isinstance(i[0], Atom):
                        if isinstance(i[1], str): # TODO we can also check the format of bond info here (what format should we use?)
                            continue
                check_pass *= 0
        if not check_pass:
            raise TypeError("Assigning wrong data type for connect. Correct data type: [[Atom(), 'bond_order_info'], ...]")

    def init_connect_in_caa(self) -> None:  # this should actually go the residue class
        """
        Initiate connectivity for this atom in a canonical amino acid. (and common solvents)
        find connect atom base on:
        1. chem.residue.RESIDUE_CONNECTIVITY_MAP
        2. parent residue name
        * Using standard Amber atom names and C/N terminal name.
          (TODO make this a standard and convert other atom name formats)
        save found list of Atom object to self._connect (make the object to a connected state)
        TODO(qz) make sure also leave space for bond order info
        """
        connect = []
        parent_residue = self.parent
        if parent_residue.name in chem.solvent.RD_SOLVENT_LIST:
            cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP[parent_residue.name][self.name]
        elif parent_residue.is_canonical():
            r = parent_residue
            r1 = parent_residue.chain[0]
            rm1 = parent_residue.chain[-1]
            if r is r1:
                # N terminal
                cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP_NTERMINAL[parent_residue.name][self.name]
            else:
                if r == rm1:
                    # C terminal
                    cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP_CTERMINAL[parent_residue.name][self.name]
                else:
                    cnt_atomnames = chem.residue.RESIDUE_CONNECTIVITY_MAP[parent_residue.name][self.name]
        else:
            _LOGGER.error(f"wrong method of getting connectivity of non-canonical residue {self.parent}. use Residue.init_connect_ncaa.")
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
                connect.append(cnt_atom)
            except ResidueDontHaveAtom as e:
                _LOGGER.warning(f"missing connecting atom {e.atom_name} of {self}. Structure maybe incomplete.")
        self._connect = connect

    #endregion

    #region === Getter-Property (ref) ===
    def radius(self, method: str = "ionic") -> float:
        """Gets the atomic radii with specified definition.
        Args:
            method: the method to determine radius the atom.
            (current available keywords)
                ionic: (ionic radius) for both metal and donor atom
                vdw: (Van Der Waals radius) for both metal and donor atom
        Returns:
            a float of the radius
        """
        radius = chem.get_atom_radii(self.element, method)
        return radius

    @dispatch
    def distance_to(self, point: Atom) -> float:
        """Get the distance to the other atom or a point."""
        return mh.get_distance(self.coord, point.coord)

    @dispatch
    def distance_to(self, point: tuple) -> float:  # pylint: disable=function-redefined
        """Get the distance to the other atom or a point."""
        if type(point) == tuple:
            return mh.get_distance(self.coord, point)
        else:
            return mh.get_distance(self.coord, point.coord)

    def attached_protons(self) -> List[Atom]:
        """find all protons attached to self"""
        result = list(filter(lambda a: a.element == "H", self.connect))
        return result

    #endregion

    #region === Check ===
    def is_donor_atom(self) -> bool:
        """check if the atom is a donor atom to a coordination center"""
        return self.name in chem.metal.DONOR_ATOM_LIST

    def is_connected(self) -> bool:
        """check if self is in the connected state"""
        return self._connect is not None

    #endregion

    #region == Special ==
    def __str__(self):
        return f"Atom({self._name}, {self._idx}, {self._coord}, {self._parent}, {self._b_factor}, {self._element}, {self._charge} )"

    #endregion

    @dispatch
    def _(self):
        """
        dummy method for dispatch
        """
        pass
