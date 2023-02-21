"""EnzyHTP is a holistic platform that allows high-throughput molecular simulation of enzymes. 
Molecular simulations, such as quantum mechanics (QM), molecular mechanics (MM), and QM/MM modeling, 
are highly applicable to the design and discovery of new biocatalysts. Molecular simulations provide 
time-resolved, atomic and electronic detail for enzymatic reactions, complementing experimental sequence 
and structure-level information. As such, expanding molecular simulation data can improve the predictive 
power of machine learning models to evaluate mutation effects in enzyme catalysis. However, large-scale 
collection of molecular simulation data presents a significant challenge due to complex demands. 
To build an enzyme model appropriate for simulations, multiple hierarchies of structural definitions and 
treatments must be established such as protein stoichiometry, binding site, predicting amino acid protonation 
state, addition of missing residues, performing an amino acid substitution, and creating reacting species. 
Most enzyme modeling practices use similar structural operations but rely on manual curation, which is 
highly inefficient and hampers reproducibility. EnzyHTP, a high-throughput enzyme simulation tool, bypasses 
these issues through automation of molecular model construction, mutation, sampling and energy calculation.
uthor: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-26
"""

from .structure import (
    Structure,
    Chain,
    Residue,
    Atom,
    Ligand,
    MetalUnit,
    Solvent,
    PDBParser,
)

from .core import (
    EnvironmentManager,
    MissingEnvironmentElement,
    InvalidResidueCode,
    _LOGGER,
    write_data,
)

from ._config import config

from .preparation import PDBLine, read_pdb_lines
from .chemical import ResidueType

from ._interface import Interface


interface = Interface(config)

from .analysis import ( electric_field )

"""Singleton interface for all softwares enzy_htp is capable of interfacing with. SHOULD NOT be used by
end users and instead exists for developers to access external software. See enzy_htp/_interface/interface.py
for full class defintion."""
