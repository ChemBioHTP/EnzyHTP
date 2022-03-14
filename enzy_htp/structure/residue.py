"""TODO DOCUMENATION"""
from __future__ import annotations

import numpy as np
from typing import Tuple
from plum import dispatch
from enzy_htp.chemical import ResidueType, METAL_MAPPER


class Residue:
    """TODO DOCUMENTATION"""
    def __init__(self, residue_key, atoms):
        self.atoms = atoms
        self.residue_key = residue_key
        (chain, name, num) = self.residue_key.split('.')
        self.chain_ = chain
        self.name = name
        self.num_ = int(num)
        self.rtype_ = ResidueType.UNKNOWN
        line_idxs = np.array(list(map(lambda aa: aa.line_idx, self.atoms)))
        self.min_line_ = np.min(line_idxs)
        self.max_line_ = np.max(line_idxs)
        # TODO what are the checks that we should be doing here?

    def __determine_residue_type( self ):
        # TODO finish this algorithm
        # 1. canoncial code => canonical
        # 2. Solvent = WAT or HOH, also solvent ions (NA+, CL-)
        # 3. Metal-center: only from metal center residue 
        #.4. Non-Canonical/Ligand => similar but ligand will be in its own chain... ligand will NEVER be in the same cahin as a canonical amino acid
        pass

    def empty_chain(self) -> bool:
        return not len(self.chain_.strip())

    def chain(self) -> str:
        return self.chain_

    def set_chain(self, val : str ) -> None:
        # TODO:CJ maybe fix this? its kinda janky
        for idx in range(len(self.atoms)):
            self.atoms[idx].chain_id = val
        self.chain_ = val

    def num(self) -> int:
        return self.num_

    def min_line(self) -> int:
        return self.min_line_
    
    def max_line(self) -> int:
        return self.max_line_

    def is_metal(self) -> bool:
        return self.name in METAL_MAPPER

    def line_range(self) -> Tuple[int,int]:
        return (self.min_line(), self.max_line())

    def neighbors(self, other : Residue) -> bool:
        """Checks if two residues from the same PDB are neighbors."""
        return (abs(self.min_line()-other.max_line()) == 1 ) or ( abs(self.max_line() - other.min_line()) == 1 )

    @dispatch
    def rtype(self, new_rtype: ResidueType) -> None:
        pass

    @dispatch
    def rtype( self ) -> ResidueType:
        pass
    #TODO add operator overloading to move the chain?
