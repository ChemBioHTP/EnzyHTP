"""Structure sub-module which encodes functionality for the StructureCluster class. This includes:

    + StructureCluster: The actual class which contains a List of geometrically similar Structure's.
    + cluster_structures: The free function which creates a List[StructureCluster] from a List[Structure].

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2025-07-06
"""

from .structure_cluster import (
    StructureCluster,
    cluster_structures 
)
