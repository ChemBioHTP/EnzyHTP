# EnzyHTP
  We built a holistic platform EnzyHTP that allows the molecular simulation of enzymes in a high-throughput manner. Molecular simulations, such as quantum mechanics (QM), molecular mechanics (MM), and QM/MM modeling, have been widely applied to guide the design and discovery of new biocatalysts. They inform the time-resolved, atomic (even electronic) detail of enzymatic reactions, which complements the sequence and structure-level information obtained in experiments. As such, augmenting molecular simulation data holds great promise to improve the predictive power of machine learning models to evaluate mutation effects for enzyme catalysis. However, large-scale collection of molecular simulation data presents a big challenge. Multiple hierarchies of structural treatments are necessary for building a simulation-ready enzyme model, including determining protein stoichiometry, identifying the binding site, predicting amino acid protonation state, adding missing residues, performing an amino acid substitution, creating reacting species, and so on. In most enzyme modeling practices, these structural operations rely on manual curation, which is highly inefficient and hampers reproducibility. Here we are developing a high-throughput enzyme simulation tool, EnzyHTP, that automates molecular model construction, mutation, sampling, and energy calculation.
![](Four_modules.png)

# Requirement
## External Program
- AmberTool/Amber
- Gaussian
- Multiwfn (for wavefunction analysis)
## Python Package
- python >= 3.6
- numpy
- pdb2pqr
- openbabel

# Installation 
## dependence
0. Install conda & create an environment
1. install numpy `conda install numpy`
2. Install openbabel `conda install openbabel -c conda-forge`
3. Install pdb2pqr 
```
git clone https://github.com/Electrostatics/pdb2pqr.git
cd pdb2pqr
pip install .
```
3. Install Multiwfn (install demo in author's blog: http://sobereva.com/454) (The LMO func seems not working for WSL)
