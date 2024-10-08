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
Complete documentation for EnzyHTP is in progress, and its code is under refactoring. 
Stay tuned for documentation and the refactored code. 
For an in-depth description of EnzyHTP, please refer to [our paper](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01424) 
and corresponding use cases under `/Test_file/FAcD_expanse`

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

## 1. Automated Installation

From the base directory of the package, i.e., the same directory as this file's.

1. Run the command `source /path/to/EnzyHTP/dev-tools/conda-install --name <YOUR_ENV_NAME>`. If you don't specify the name, it will be installed to your current conda environment.
2. (optional) Install Gaussian16.
3. (optional) Install Multiwfn (install demo in author's [blog](http://bbs.keinsci.com/thread-12020-1-1.html): (The LMO func seems not working for WSL) (Note that run Multiwfn on ACCRE requires loading the GCC module).

## 2. Automated Installation as a Dependency.

As we are to install the EnzyHTP as a package, the installation of Gaussian16 and Multiwfn are omitted.

Run the command `source /path/to/EnzyHTP/dev-tools/conda-install --name <YOUR_ENV_NAME>`. If you don't specify the name, it will be installed to your current conda environment.

## 3. Manual Installation

1. Install conda & (optional) create an environment
2. install numpy `conda install numpy`
3. Install openbabel `conda install openbabel -c conda-forge`
4. Install pdb2pqr 
```
git clone https://github.com/Electrostatics/pdb2pqr.git
cd pdb2pqr
pip install .
```
5. Install [AmberTool](https://ambermd.org/AmberTools.php)
6. (optional) Install RDKit `conda install -c conda-forge rdkit`
7. (optional) Install Gaussian16.
8. (optional) Install Multiwfn (install demo in author's [blog](http://bbs.keinsci.com/thread-12020-1-1.html): (The LMO func seems not working for WSL) (Note that run Multiwfn on ACCRE requires loading the GCC module).
