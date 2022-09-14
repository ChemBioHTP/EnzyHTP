'''Definition for the Atom class. Meant to be the base unit of structural information: coordinate and
topology (atom type & connectivity). Serve solely for storing and accessing data. 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
'''
from __future__ import annotations
from typing import Union, Any

import numpy as np
import pandas as pd


class Atom:
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
        self._parent = parent # could to be add later
        # optional
        self._b_factor = ds['b_factor']
        self._element = ds['element_symbol'].strip()
        if not np.isnan(ds['charge']):
            self._charge = float(ds['charge'])
