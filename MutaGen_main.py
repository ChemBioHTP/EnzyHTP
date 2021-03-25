from Class_PDB import *
from Class_Conf import *

#This is the main file
def main():
    #initial PDB
    #use tleap to randomly mutate the PDB
    for i in range(1):
        PDB1=PDB('2kz2init_amb.pdb',wk_dir='./2kz2')
        PDB1.Add_MutaFlag('XA12G')
        PDB1.PDB2PDBwLeap()
    #use minimization to relax each mutated PDB
        PDB1.PDB2FF()
        PDB1.PDBMin()
    #run MD
        PDB1.rm_wat()
        PDB1.PDB2FF()

        PDB1.conf_prod['nstlim'] = 100000000 # Edit MD configuration (see default in Class_Conf.py - Config.Amber)
        PDB1.PDBMD(tag='XA12G')

        #print(Check_PDB(PDB3_path,Flag), file='min_cache/Check'+Flag+'.log')

if __name__ == "__main__":
  main()


