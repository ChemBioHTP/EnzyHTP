# EnzyHTP
we built a holistic platform EnzyHTP that allows the collection of both dynamics from MD and electronic structure data from QM or QM/MM in a high-throughput manner.

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
