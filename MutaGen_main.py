from Class_PDB import *
from Class_Conf import *

#This is the main file

#initial PDB
#use tleap to randomly mutate the PDB
for i in range(2):
    PDB1=PDB('2kz2init_amb.pdb')
    PDB1.Random_MutaFlag()
    PDB1.PDB2PDBwLeap()
#use minimization to relax each mutated PDB
    PDB1.PDB2FF()
    PDB1.PDBMin()
#run MD
    os.system('mkdir MD')
    PDB1.rm_wat()
    PDB1.PDB2FF()
    MD=Conf()
    MD.set_min()
    MD.set_heat()
    MD.set_equi()
    MD.set_prod()
    MD.deploy(path='./MD')
    PDB1.PDBMD(MD.deploy_path)

    #print(Check_PDB(PDB3_path,Flag), file='min_cache/Check'+Flag+'.log')

