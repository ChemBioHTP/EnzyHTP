'''Definition for the Atom class. Meant to be the base unit of structural information: coordinate and
topology (atom type & connectivity). Serve solely for storing and accessing data. 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
'''
from __future__ import annotations
from typing import Union, Any

import numpy as np


class Atom:
    '''Base unit of structure in enzy_htp. Store coordinate and topology information
    in a Structure() object. Serve solely for storing and accessing data.
    * not specific to any file format that store a structure

    Attributes:
        (nessessary)
        name : The name of the atom as a string. often refer to a specific connectivity.
        index: An integer number for the atomic index.
        coord: (x,y,z) for cartesian coordinate of the atom.
        parent/residue: the parent residue that this atom belongs to.
        (optional)
        b_factor : A float representing the temperature factor (b factor).
        charge : Charge of the atom.
        element : Character representing element.
    '''

    def __init__(self, **kwargs): #@shaoqz: @imp we should have a seperate IO module parsing PDB and have functions there as constructor and exportor of a structure object.
        '''Constructor meant to be called as Atom(**row), where row is of type pandas.Series'''

        for k, v in kwargs.items():
            setattr(self, k, v)

        if not np.isnan(self.charge):
            self.charge = int(self.charge) #@shaoqz: charge should be float.