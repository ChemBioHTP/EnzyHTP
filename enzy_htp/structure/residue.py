"""TODO DOCUMENATION"""

from plum import dispatch
from enzy_htp.chemical import ResidueType


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
        # TODO what are the checks that we should be doing here?

    def __determine_residue_type( self ):
        # TODO finish this algorithm
        # 1. canoncial code => canonical
        # 2. Solvent = WAT or HOH, also solvent ions (NA+, CL-)
        # 3. Metal-center: only from metal center residue 
        #.4. Non-Canonical/Ligand => similar but ligand will be in its own chain... ligand will NEVER be in the same cahin as a canonical amino acid
        pass

    def chain(self):
        return self.chain_

    def num(self):
        return self.num_
	
    @dispatch
    def rtype(self, new_rtype: ResidueType) -> None:
        pass
	
    @dispatch
    def rtype( self ) -> ResidueType:
        pass
    #TODO add operator overloading to move the chain?
