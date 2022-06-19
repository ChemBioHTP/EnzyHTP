"""
TODO(CJ)
"""
from typing import Dict
from copy import deepcopy

import enzy_htp.structure as es

class Restrictions:
    """
    """
    def __init__(self, mapper):
        self.mapper = deepcopy(mapper) 



    def update(self):
        pass


    def lock_chain(self, cname:str) -> None:
        pass

    def lock_residue(self, cname:str) -> None:
        pass

    def lock_residues(self, cname:str) -> None:
        pass

 
    def apply(self, 


def default_restriction_dict() -> Dict:
    """
    """
    result = dict()
	result['locked'] = False
	result['illegal_targets'] = deepcopy([])
    result['no_size_increase'] = False
    result['no_size_decrease'] = False
	result['no_polarity_change'] = False
	result['force_polarity_change'] = False
    return deepcopy(result)

def restriction_object(pdb:str) -> Restrictions:
    """TODO(CJ)"""
	struct: es.Structure = es.structure_from_pdb(pdb)
    pass
