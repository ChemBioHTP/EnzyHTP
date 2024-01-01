"""Definition for the Atom class. Meant to be the base unit of structural information: coordinate and
topology (atom type & connectivity). Serve solely for storing and accessing data.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
from __future__ import annotations
import re
import sys
from typing import Any, List, Tuple, Union, Dict
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
        atom_type: a non-standard attribute of Atom(). It is specific to a force field or software.
        connect : the connectivity of this atom. (a list of reference of connected Atom objs)
    """

    def __init__(self, ds: Union[Dict, pd.Series], parent=None):
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
        self._atom_type = None
        ds_keys = ds.keys()
        if "atom_number" in ds_keys and not np.isnan(ds["atom_number"]):
            self._idx = ds["atom_number"]
        if "b_factor" in ds_keys and not np.isnan(ds["b_factor"]):
            self._b_factor = ds["b_factor"]
        if "element_symbol" in ds_keys and ds["element_symbol"].strip() != "":
            self._element = ds["element_symbol"].strip()
        if "charge" in ds_keys and not np.isnan(ds["charge"]):
            self._charge = float(ds["charge"])
        if "atom_type" in ds_keys and ds["atom_type"].strip() != "": #TODO(qz): find a way to unify
            self._atom_type = ds["atom_type"]

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
            elif self.parent is not None and self.parent.is_metal():
                return self.parent.element
            else:
                # case: in ligand atoms are named like this H1/1H
                # case: CH3 from FAH
                # case: CL1 from DCE TODO (add a mapper for all 2-letter elements, 
                # and give warning if it can be proven that current info is not enough
                # to deduce the ele. e.g.: is CO1 C or Co)
                clean_name = self.name.lstrip('0123456789').rstrip('0123456789')
                return re.match("^[A-Z][a-z]?", clean_name).group()           

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
    def connect(self) -> List[Tuple[Atom, str]]:
        """getter for _connect, the list for (atoms, bond_type) it connects"""
        if not self.is_connected():
            _LOGGER.error(f"There are no connection info for {self}. "
                            "Please initiate it use structure.structure_operation.init_connectivity()")
            raise AttributeError
        return self._connect

    @connect.setter
    def connect(self, val: List[Tuple[Atom, str]]):
        """setter of connect is nessessary for residue level objects to set connectivity for atoms
        raise an error if the assigned val doesn't align with the format"""        
        # check for data format
        type(self)._check_assigning_connect_data_type(val)

        # assign the value if passed
        self._connect = val

    def connect_to(self, other: Atom, bond: str=None):
        """connect {self} with {other}. add each other to .connect with {bond}"""
        # connect self to other
        if self._connect is None:
            self._connect = []

        if other in self.connect_atoms:
            _LOGGER.debug(f"{other} already in connect of {self}. Doing nothing.")
            return

        connect_info = (other, bond)
        type(self)._check_connect_element_data_type(connect_info)
        self.connect.append(connect_info)

        # connect other to self
        if other._connect is None:
            other._connect = []

        if self in other.connect_atoms:
            _LOGGER.debug(f"{self} already in connect of {other}. Doing nothing.")
            return

        connect_info = (self, bond)
        type(other)._check_connect_element_data_type(connect_info)
        other.connect.append(connect_info)

    @classmethod
    def _check_assigning_connect_data_type(cls, val: Any) -> bool:
        """Function used by the setter of connect. Check if the data format of val meetings the definiation.
        The format is defined in this function.
        Raise an error if the format is not met.
        Format: [(Atom(), "bond_order_info"), ...]"""
        check_pass = 1
        if isinstance(val, list):
            for i in val:
                if cls._check_connect_element_data_type(i, raise_error=False):
                    continue
                check_pass *= 0
        if not check_pass:
            raise TypeError("Assigning wrong data type for connect. Correct data type: [[Atom(), 'bond_order_info'], ...]")

    @classmethod
    def _check_connect_element_data_type(cls, val: Any, raise_error: bool=True) -> bool:
        """Function used by the setter of connect. Check if the data format of val meetings the definiation.
        The format is defined in this function.
        Raise an error if the format is not met.
        Format: (Atom(), "bond_order_info")"""
        check_pass = False
        if isinstance(val, tuple):
            if isinstance(val[0], Atom):
                if isinstance(val[1], str) or (val[1] is None):
                    check_pass = True
        if raise_error:
            if not check_pass:
                raise TypeError("Adding wrong data type to connect. Correct data type: (Atom(), 'bond_order_info')")
        else:
            return check_pass

    @property
    def key(self) -> str: # TODO change the name to key_str
        """Gets the Atom()'s key which can be used in conjuection with the Structure.get_atom method.
        Format is <chain_name>.<residue_index>.<atom_name>. If the atom does not have a parent residue or
        chain, an empty value is given (e.g. "..CA")
        """
        
        tokens:List[str] = ["", "", self.name]

        if self.parent is not None:
            tokens[1] = str(self.parent.idx)
            
            if self.parent.parent is not None:
                tokens[0] = self.parent.parent.name
    
        return '.'.join(tokens)
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
    def distance_to(self, point: Tuple) -> float:  # pylint: disable=function-redefined
        """Get the distance to the other atom or a point."""
        return mh.get_distance(self.coord, point)
    
    @dispatch
    def angle_with(self, point_1: Atom, point_2: Atom) -> float:
        """Get the angle to the other 2 atoms or 2 points."""
        return mh.get_angle(self.coord, point_1.coord, point_2.coord)

    @dispatch
    def angle_with(self, point_1: Tuple, point_2: Tuple) -> float: # pylint: disable=function-redefined
        """Get the angle to the other 2 atoms or 2 points."""
        return mh.get_angle(self.coord, point_1, point_2)

    @dispatch
    def dihedral_with(self, point_1: Atom, point_2: Atom, point_3: Atom) -> float:
        """Get the dihedral to the other 3 atoms or 3 points."""
        return mh.get_dihedral(self.coord, point_1.coord, point_2.coord, point_3.coord)

    @dispatch
    def dihedral_with(self, point_1: Tuple, point_2: Tuple, point_3: Tuple) -> float: # pylint: disable=function-redefined
        """Get the dihedral to the other 3 atoms or 3 points."""
        return mh.get_dihedral(self.coord, point_1, point_2, point_3)

    def attached_protons(self) -> List[Atom]:
        """find all protons attached to self"""
        result = list(filter(lambda a: a[0].element == "H", self.connect))
        result = [cnt_atom for cnt_atom, bond in result]
        return result

    @property
    def connect_atoms(self) -> List[Atom]:
        """get all atoms of connection without the bond info"""
        if not self.is_connected():
            _LOGGER.error(f"There are no connection info for {self}. "
                            "Please initiate it use structure.structure_operation.init_connectivity()")
            raise AttributeError
        return [atom for atom, bond in self._connect]
    
    @property
    def connect_atom_names(self) -> List[str]:
        """get all atom names of connection."""
        return [atom.name for atom in self.connect_atoms]

    #endregion

    #region === Check ===
    def is_donor_atom(self) -> bool:
        """check if the atom is a donor atom to a coordination center"""
        return self.name in chem.metal.DONOR_ATOM_LIST

    def is_connected(self) -> bool:
        """check if self is in the connected state"""
        return self._connect is not None
    
    def is_mainchain_atom(self) -> bool:
        """check if self is a mainchain atom.
        (i.e.: N CA C)"""
        if self.residue.is_canonical():
            if self.name in "N CA C".split():
                return True
        elif self.residue.is_modified():
            raise Exception("TODO prob determine start/end atom and deduce mainchain")

        return False

    #endregion

    #region == Special ==
    def __str__(self):
        return f"Atom({self._name}, {self._idx}, {self._coord}, {self._parent}, {self._b_factor}, {self._element}, {self._charge} )"

    def __repr__(self):
        return f"<{self} at {hex(id(self))}>"
    #endregion

    @dispatch
    def _(self):
        """
        dummy method for dispatch
        """
        pass
