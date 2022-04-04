"""Specialization of the Residue() class for a Ligand. 

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-04-03
"""

import numpy as np

from .atom import Atom
from typing import List
from .residue import Residue
from enzy_htp.chemical import enum as renum

from enzy_htp.core import file_system as fs
from enzy_htp.core import (
    UnsupportedFileType,
    _LOGGER,
)


class Ligand(Residue):
    """Represents a specific Ligand found in a .pdb file.

    Attributes:
		self.net_charge : The net charge of the molecule.
    """
    def __init__(self, residue_key: str, atoms: List[Atom], **kwargs):
        """Constructor for Ligand. Identical to Residue() ctor but also takes net_charge value."""
        self.net_charge = kwargs.get("net_charge", None)
        Residue.__init__(self, residue_key, atoms)
        self.set_rtype( renum.ResidueType.LIGAND )

    def is_ligand(self) -> bool:
        """Checks if the Residue is a ligand. Always returns True for this specialization."""
        return True

    def set_residue_number(self, num: int) -> None:
        """Changes the resdiue number for all of the substituent atoms."""
        for idx, aa in enumerate(self.atoms):
            self.atoms[idx].residue_number = num

    def get_net_charge(self) -> int:
        """Getter for the net_charge attribute."""
        return self.net_charge

    def build(self, out_path: str) -> None:
        """Method that builds the given ligand to the specified path, making sure it is a pdb filepath."""
        ext = fs.get_file_ext(out_path).lower()
        if ext != ".pdb":
            _LOGGER.error(f"The supplied file path '{out_path}' does not ahve a '.pdb' extension. Exiting...")
            exit( 1 )
        lines = list(
            map(
                lambda pr: pr[1].to_pdb_line(a_id=pr[0] + 1, c_id=" "),
                enumerate(self.atoms),
            )
        ) + ["TER", "END"]
        fs.write_lines(out_path, lines)

def residue_to_ligand(ptr: Residue, net_charge : float = None) -> Ligand:
    """Convenience function that converts Residue to ligand."""
    return Ligand(ptr.residue_key, ptr.atoms, net_charge=net_charge)
