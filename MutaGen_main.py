from LeapPDB import *
from AmberMaps import *
from TestTools import *

#This is the main file

#initial PDB
PDB1_path='2kz2init_amb.pdb'
#use tleap to randomly mutate the PDB
for i in range(100):
    Flag=FlagGen(PDB1_path)
    PDB2_path=PDB2Leap(PDB1_path, Flag)
#use minimization to relax each mutated PDB
    PDB3_path=PDBMin(PDB2_path)

