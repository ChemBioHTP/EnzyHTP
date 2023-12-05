





from typing import List, Tuple, Union


import enzy_htp.interface as interface

from ..structure import Structure, Residue



class StructureRegion(Structure):
    pass

    



def create_structure_resion(stru:Structure,
                            selection:List[Union[Tuple[str,int],int]],
                            nterm_cap:str=None,
                            cterm_cap:str=None
                            ) -> StructureRegion:

 
    
    pass
