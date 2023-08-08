import datetime
import sys
import os
import glob
import numpy as np
from Class_PDB import *
from Class_Conf import *
from Class_ONIOM_Frame import *
from helper import write_data, line_feed

# settings
Config.n_cores = 8
Config.max_core = 2000
#Config.PC_cmd = 'srun' # srun does not work like this use mpi instead
Config.Amber.conf_equi['nstlim'] = 500000
Config.Amber.conf_prod['nstlim'] = 55000000
Config.Amber.conf_prod['ntwx'] = '50000'
Config.debug = 1
wkflow_log_path = './Gen.py.log'

def main():
        starttime = datetime.datetime.now()
        of = open(wkflow_log_path, 'wb', buffering=0)
        # deploy mutation
        of.write(('Working on: '+ 'C131V' +line_feed).encode('utf-8'))
        # --- Preparation ---
        pdb_obj = PDB('/scratch/jiany37/HMT/ALI/input_prep/HMT_ALI_RDock_aH_amber_new.pdb')
        of.write(('Preparation: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))
        pdb_obj.conf_min['DISANG']=pdb_obj.dir+'/MD/0.rs'
        pdb_obj.conf_heat['DISANG']=pdb_obj.dir+'/MD/0.rs'
        pdb_obj.conf_equi['DISANG']=pdb_obj.dir+'/MD/0.rs'
        pdb_obj.conf_prod['DISANG']=pdb_obj.dir+'/MD/0.rs'

        # --- Operation ---
        # Mutation
        pdb_obj.Add_MutaFlag('C131V')
        pdb_obj.PDB2PDBwLeap()
        of.write(('Mutation: p2pwl: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))
        pdb_obj.rm_allH()
        pdb_obj.get_protonation(if_prt_ligand=0)

        # use minimization to relax each mutated PDB
        pdb_obj.PDB2FF(local_lig=0)
        pdb_obj.PDBMin(engine='Amber_GPU')
        pdb_obj.rm_wat()
        of.write(('Mutation: PDBMin: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))
        # -- DISANG setup --
        pdb_obj.conf_min['nmropt_rest'] = '1'
        pdb_obj.conf_heat['nmropt_rest'] = '1'
        pdb_obj.conf_equi['nmropt_rest'] = '1'
        pdb_obj.conf_prod['nmropt_rest'] = '1'
        pdb_obj.get_stru()
        a1=str(pdb_obj.stru.ligands[0].S1.id)
        a2=str(pdb_obj.stru.ligands[1].C1.id)
        a3=str(pdb_obj.stru.ligands[1].I1.id)
        pdb_obj.conf_prod['rs_constraints']=[
                               {
                               'iat':[a1,a2],'r1': '1.00','r2': '2.00','r3': '3.80','r4': '5.00',
                               'rk2':'0.0','rk3':'100.0','ir6':'1','ialtd':'0'
                                },
                               {
                               'iat':[a1,a2,a3],'r1': '90.0','r2': '150.0','r3': '210.0','r4': '270.0',
                               'rk2':'100.0','rk3':'100.0','ir6':'1','ialtd':'0'
                                },
                            ]
        # --- Sample with MD ---
        pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1)
        pdb_obj.PDBMD(engine='Amber_GPU')
        of.write(('MD: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))
        # --- QM cluster ---
        pdb_obj.nc2mdcrd(start=101,step=10)
        atom_mask = ':217,218'
        g_route = '# pbe1pbe/def2TZVP nosymm'
        pdb_obj.PDB2QMCluster(atom_mask, g_route=g_route, ifchk=0)
        of.write(('QM Cluster: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))

        endtime = datetime.datetime.now()
        print(endtime - starttime)
        of.close()

if __name__ == "__main__":
        main()