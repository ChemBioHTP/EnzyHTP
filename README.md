# enzyme_workflow
Using AlleyCat as a model enzyme, we will develop Python 3 libraries that facilitate the enzyme mutant model construction, input file creation, job monitor, output data collection, and analysis.

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
1. Install openbabel `conda install openbabel -c conda-forge`
2. Install pdb2pqr 
```
git clone https://github.com/Electrostatics/pdb2pqr.git
cd pdb2pqr
pip install .
```
3. Install Multiwfn (install demo in author's blog: http://sobereva.com/454)
