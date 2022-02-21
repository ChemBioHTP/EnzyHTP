import numpy as np


class Atom:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        # self.name = kwargs.get('atom_name',str())
        ##self.x = kwargs.get('x_coord',float())
        ##self.y = kwargs.get('y_coord',float())
        ##self.z = kwargs.get('z_coord',float())
        # self.atom_number = kwargs.get('atom_number', int())
        # self.residue_number = kwargs.get('residue_number', int())
        # self.residue_name = kwargs.get('residue_name', int())
        # self.chain_id = kwargs.get('chain_id', str())

    def __add__(self, offset):
        return Atom(
            name=self.name,
            x=self.x_coord + other,
            y=self.y_coord + other,
            z=self.z_coord + other,
        )

    def __sub__(self, offset):
        return Atom(
            x_coord=self.x - other,
            y_coord=self.y - other,
            z_coord=self.z - other,
            **self.__dict__,
        )

    def residue_key(self):
        return f"{self.chain_id}.{self.residue_name}.{self.residue_number}"

    def distance(self, other):
        return np.sqrt(
            (self.x_coord - other.x_coord) ** 2
            + (self.y_coord - other.y_coord) ** 2
            + (self.z_coord - other.z_coord) ** 2
        )

    def __getitem__(self, key):
        return self.__dict__[key]
