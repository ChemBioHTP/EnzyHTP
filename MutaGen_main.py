from Class_PDB import *
from Class_Conf import *

#This is the main file
def main():
    #initial PDB
    #use tleap to randomly mutate the PDB
    for i in range(1):
        PDB1=PDB('2kz2init_amb.pdb')
        PDB1.Add_MutaFlag('random')
        PDB1.PDB2PDBwLeap()
    #use minimization to relax each mutated PDB
        PDB1.PDB2FF()
        PDB1.PDBMin()
    #run MD
        os.system('mkdir MD')
        PDB1.rm_wat()
        PDB1.PDB2FF()

        #Maybe these should go into the PDBMD func?
        # change the structure of conf to ade-like
        MD=Conf()
        MD.set_MD_min()
        MD.set_MD_heat()
        MD.set_MD_equi()
        MD.set_MD_prod()
        MD.deploy(path='./MD')
        PDB1.PDBMD(conf_path=MD.deploy_path)

        #print(Check_PDB(PDB3_path,Flag), file='min_cache/Check'+Flag+'.log')

if __name__ == "__main__":
  main()


