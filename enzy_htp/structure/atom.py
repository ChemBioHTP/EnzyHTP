"""Definition for the Atom class. Meant to be the base unit of structural information: coordinate and
topology (atom type & connectivity). Serve solely for storing and accessing data. 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
from __future__ import annotations
import sys
from typing import Tuple

import numpy as np
import pandas as pd
from enzy_htp.core.doubly_linked_tree import DoubleLinkNode
from enzy_htp.core.logger import _LOGGER


class Atom(DoubleLinkNode):
    """Base unit of structure in enzy_htp. Store coordinate and topology information
    in a Structure() object. Serve solely for storing and accessing data.
    * not specific to any file format that store a structure

    Attributes:
        (nessessary)
        name : The name of the atom as a string. often refer to a specific connectivity.
        coord: (x,y,z) for cartesian coordinate of the atom.
        parent/residue: the parent residue that this atom belongs to.
        (optional)
        idx: An integer number for the atomic index.
        b_factor : A float representing the temperature factor (b factor).
        charge : Charge of the atom.
        element : Character representing element.
    """

    def __init__(self, ds: pd.Series, parent = None):
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
        ds_keys = ds.keys()
        if "atom_number" in ds_keys:
            self._idx = ds["atom_number"]
        if "b_factor" in ds_keys:
            self._b_factor = ds["b_factor"]
        if "element_symbol" in ds_keys:
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
        if self._idx == None:
            _LOGGER.warning(f"{sys._getframe().f_back.f_back.f_code.co_name}:{sys._getframe().f_back.f_code.co_name} is using index of an non-indexing atom: {self}")
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
        """synonym for element"""
        return self._element
    @element.setter
    def element(self, val):
        self._element = val

    @property
    def charge(self):
        """synonym for charge"""
        return self._charge
    @charge.setter
    def charge(self, val):
        self._charge = val
    #endregion

    #region == Special ==
    def __str__(self):
        return f"Atom({self._name}, {self._idx}, {self._coord}, {self._parent}, {self._b_factor}, {self._element}, {self._charge} )"
    #endregion

