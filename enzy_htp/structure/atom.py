'''Definition for the Atom class. Meant to be the base unit of structural information: coordinate and
topology (atom type & connectivity). Serve solely for storing and accessing data. 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
'''
from __future__ import annotations
from typing import Tuple

import numpy as np
import pandas as pd
from enzy_htp.core.doubly_linked_tree import DoubleLinkNode


class Atom(DoubleLinkNode):
    '''Base unit of structure in enzy_htp. Store coordinate and topology information
    in a Structure() object. Serve solely for storing and accessing data.
    * not specific to any file format that store a structure

    Attributes:
        (nessessary)
        name : The name of the atom as a string. often refer to a specific connectivity.
        idx: An integer number for the atomic index.
        coord: (x,y,z) for cartesian coordinate of the atom.
        parent/residue: the parent residue that this atom belongs to.
        (optional)
        b_factor : A float representing the temperature factor (b factor).
        charge : Charge of the atom.
        element : Character representing element.
    '''

    def __init__(self, ds: pd.Series, parent = None):
        '''Constructor of Atom(), where ds is of type pandas.Series'''
        # nessessary
        self._name = ds['atom_name'].strip()
        self._idx = ds['atom_number']
        self._coord = (ds['x_coord'], ds['y_coord'], ds['z_coord'])
        self.set_parent(parent)
        self.set_ghost_children()
        # optional
        self._b_factor = ds['b_factor']
        self._element = ds['element_symbol'].strip()
        if not np.isnan(ds['charge']):
            self._charge = float(ds['charge'])

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
        """synonym for parent"""
        return self._b_factor
    @b_factor.setter
    def b_factor(self, val):
        self._b_factor = val

    @property
    def element(self):
        """synonym for parent"""
        return self._element
    @element.setter
    def element(self, val):
        self._element = val

    @property
    def charge(self):
        """synonym for parent"""
        return self._charge
    @charge.setter
    def charge(self, val):
        self._charge = val
    #endregion
