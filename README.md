[![DOI](https://zenodo.org/badge/459668860.svg)](https://zenodo.org/badge/latestdoi/459668860)
[![build](https://github.com/google/yapf/actions/workflows/ci.yml/badge.svg)](https://github.com/google/yapf/actions)

# EnzyHTP
EnzyHTP is a holistic platform that allows high-throughput molecular simulation of enzymes. 
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
![](resources/four_modules_whitebg.png)

# Documentation
The documents are in the making at https://enzyhtp-doc.readthedocs.io/en/latest/
For an in-depth description of EnzyHTP, you can also refer to [our paper](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01424) 
and corresponding use cases (tests) under `/test` and `/test_integration`

# Installation

https://enzyhtp-doc.readthedocs.io/en/latest/install.html#

# License

For non-commercial use, please find our “Non-Exclusive Non-Commercial Academic Software License Agreement” in `license.md`.

As a summary, for non-commercial use, you are allowed to (for free)
- copy
- modify

you are not allowed to:
- redistribute a modified version


For commercial use, please send a request to Vanderbilt (cttc@vanderbilt.edu) and cc Dr. Yang (zhongyue.yang@vanderbilt.edu) and QZ (qz.shao@outlook.com) for a commerical license, it will be similar to the EULA but distributed through the CTTC office.

As a summary, for commercial use, after obtaining the license from Vanderbilt, you are allowed to (for free)
- copy
- modify

you are not allowed to:
- redistribute a modified version
