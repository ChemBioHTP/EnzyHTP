import datetime
import enum
import os
from Class_PDB import *
from Class_ONIOM_Frame import *

def main():
    starttime = datetime.datetime.now()

    pdb_obj = PDB('2kz2init_amb.pdb')
    pdb_obj.PDB2FF(ifsavepdb=1)            
    pdb_obj.PDBMD(engine='Amber_pmemd_gpu') # MD
    pdb_obj.nc2mdcrd(point=100)             # sample
    # make ONIOM template
    pdb_obj.layer_atoms=['1-9','10-L']
    gout_path = pdb_obj.PDB2QMMM(qm='g16', ifchk=0)
    
    # analyze with the ensemble


    endtime = datetime.datetime.now()
    print(endtime - starttime)

if __name__ == "__main__":
    main()
