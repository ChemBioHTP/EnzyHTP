from Class_PDB import *
from Class_Conf import *

#This is the main file
def main():
    #initial PDB
    #use tleap to randomly mutate the PDB
    for i in range(1):
        PDB1=PDB('2kz2init_amb.pdb',wk_dir='./2kz2')
        print(PDB1.Add_MutaFlag('r'))
        PDB1.PDB2PDBwLeap()
    # remove causal protonation from leap (Not sure if it's better to be before Min or not)
        PDB1.rm_allH()
        PDB1.get_protonation()
    #use minimization to relax each mutated PDB
        PDB1.PDB2FF()
        PDB1.PDBMin()
    #run MD
        PDB1.rm_wat()
        PDB1.PDB2FF()
        PDB1.conf_prod['nstlim'] = 100000000 # Edit MD configuration (see default in Class_Conf.py - Config.Amber)
        PDB1.PDBMD(tag='XA12G')
    
    # demo: reset_MD_conf
        label = PDB1.Add_MutaFlag('r')
        PDB1.PDB2PDBwLeap()
        PDB1.PDB2FF(renew_lig=1) # not nessessary but available
        PDB1.PDBMin()
        PDB1.rm_wat()
        PDB1.PDB2FF()
        PDB1.reset_MD_conf()
        PDB1.PDBMD(tag=label)

        #print(Check_PDB(PDB3_path,Flag), file='min_cache/Check'+Flag+'.log')

if __name__ == "__main__":
  main()


