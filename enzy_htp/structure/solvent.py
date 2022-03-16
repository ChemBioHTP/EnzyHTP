from .atom import Atom
from .residue import Residue
from typing import List


class Solvent(Residue):
    """
    -------------
    initilize from
    PDB:        Residue.fromPDB(atom_input, input_type='PDB_line' or 'line_str' or 'file' or 'path')
    raw data:   Residue(atom_name, coord)
    Residue:    Ligand.fromResidue(atom_obj)
    -------------
    id
    name
    atoms = [atom, ...]
    parent # the whole stru
    -------------
    Method
    -------------
    set_parent
    -------------
    """

    def __init__(self, residue_key: str, atoms: List[Atom]):
        Residue.__init__(residue_key, atoms)

    def sort(self):
        pass


def residue_to_solvent(ptr: Residue) -> Solvent:
    # TODO documentation
    return Solvent(residue_key, atoms)
