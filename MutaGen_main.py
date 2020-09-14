from LeapPDB import *
from AmberMaps import *
from TestTools import *

#This is the main file

#initial PDB
PDB1_path='2kz2init_amb.pdb'
#use tleap to mutate the PDB
PDB2_path=PDB2Leap(PDB1_path, 'E37K')
#use minimization to relax the mutated PDB
PDB3_path=PDBMin(PDB2_path)

