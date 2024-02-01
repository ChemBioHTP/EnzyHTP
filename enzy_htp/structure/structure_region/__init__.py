"""The StructureRegion() class is making a sub-section of Structure(), usually to do a QM calculation, either
single point or geometry optimization, or multiscale modelling with QM/MM methods. The StructureRegion() also 
has the ability to cap the backbones of amino acids with various capping groups, managed by ResidueCap() classes.
The StructureRegion() does not own the atoms involved, and ownership is maintained by the original Structure(). 

except for capping atoms. 


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-01-22
"""


from .api import (
    StructureRegion,
    ResidueCap,
    create_region_from_selection_pattern,
    create_region_from_residue_keys,
    create_region_from_full_stru,
)
