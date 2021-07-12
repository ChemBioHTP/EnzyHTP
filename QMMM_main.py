import datetime
import enum
import os
from Class_PDB import *
from Class_ONIOM_Frame import *

def main():
    starttime = datetime.datetime.now()

    a = PDB('2kz2init_amb.pdb')
    a.PDB2FF(ifsavepdb=1)
    nc_path = a.PDBMD()
    a.layer_atoms=['1-9','10-L']
    g_temp_path = a.PDB2QMMM()

    #convert nc to mdcrd (pytraj)
    frames = Frame.fromMDCrd(nc_path)
    gjf_paths = []
    for i, frame in enumerate(frames):
        gjf_paths.append(frame.write_to_template(g_temp_path, index = str(i)))
        
    # run/qsub Gaussian job
    out_paths = Run_QMMM(gjf_paths)
    
    # analyze with the ensemble


    endtime = datetime.datetime.now()
    print(endtime - starttime)

if __name__ == "__main__":
    main()
