from core.clusters.accre import Accre
import datetime
from Class_PDB import PDB
from Class_Conf import Config
from Class_ONIOM_Frame import *
from helper import write_data

# settings
Config.n_cores = 24 # still submit this script to 24 cores used by Multiwfn
# run(f"sed -i 's/nthreads= *[0-9][0-9]*/nthreads=  {Config.n_cores}/' {Config.Multiwfn.DIR}/settings.ini", check=True, text=True, shell=True, capture_output=True)
Config.max_core = 2000
Config.PC_cmd = 'srun'
Config.Amber.conf_equi['nstlim'] = 500000
Config.Amber.conf_prod['nstlim'] = 55000000
Config.Amber.conf_prod['ntwx'] = '50000'
Config.debug = 1
data_output_path = './Mutation.dat'

def main():
    starttime = datetime.datetime.now()
    # --- Preparation ---
    pdb_obj = PDB('/scratch/jiany37/HMT/ALI/input_prep/HMT_ALI_RDock_aH_amber_new.pdb')
    
    pdb_obj.conf_min['DISANG']=pdb_obj.dir+'/MD/0.rs'
    pdb_obj.conf_heat['DISANG']=pdb_obj.dir+'/MD/0.rs'
    pdb_obj.conf_equi['DISANG']=pdb_obj.dir+'/MD/0.rs'
    pdb_obj.conf_prod['DISANG']=pdb_obj.dir+'/MD/0.rs'

    # --- Operation ---
    # Mutation
    pdb_obj.Add_MutaFlag('C131V')
    pdb_obj.PDB2PDBwLeap()
    pdb_obj.rm_allH()
    pdb_obj.get_protonation(if_prt_ligand=0)
    # use minimization to relax each mutated PDB
    pdb_obj.PDB2FF(local_lig=0)
    pdb_obj.PDBMin( engine='Amber_GPU', 
                    if_cluster_job=1,
                    cluster=Accre(),
                    period=10,
                    res_setting={'partition': 'maxwell', 
                                 'account':'csb_gpu_acc'} ) # see full setting in Conf
    pdb_obj.rm_wat()

    endtime_1 = datetime.datetime.now()
    print(f'mutation: {endtime_1 - starttime}')
    # --- Sample with MD ---
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
    pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1)
    pdb_obj.PDBMD(  engine='Amber_GPU', 
                    equi_cpu=0, 
                    if_cluster_job=1,
                    cluster=Accre(),
                    period=30,
                    res_setting={'account':'csb_gpu_acc'} )# you may also change the partition by add a 'partition' : 'xxx' to this dict (default: pascal)
    endtime_2 = datetime.datetime.now()
    print(f'MD: {endtime_2 - endtime_1}')
    # --- QM cluster ---
    pdb_obj.nc2mdcrd(start=101,step=10)
    atom_mask = ':217,218'
    g_route = '# pbe1pbe/def2TZVP nosymm'
    pdb_obj.PDB2QMCluster(  atom_mask, 
                            g_route=g_route,
                            ifchk=1,
                            if_cluster_job=1, 
                            cluster=Accre(), 
                            job_array_size=20, # this is how many job can on queue the same time
                            period=30,
                            res_setting={'account':'yang_lab'} )
    pdb_obj.get_fchk(keep_chk=0)

    endtime_3 = datetime.datetime.now()
    print(f'QMCluster: {endtime_3 - endtime_2}')
    # --- Analysis ---
    # targeting C-I bond
    a1 = int(pdb_obj.stru.ligands[1].C1)
    a2 = int(pdb_obj.stru.ligands[1].I1)
    a1qm = pdb_obj.qm_cluster_map[str(a1)]
    a2qm = pdb_obj.qm_cluster_map[str(a2)]
    # Field Strength (MM)
    E_atom_mask = ':1-216'
    Es = pdb_obj.get_field_strength(E_atom_mask, a1=a1 ,a2=a2 ,bond_p1='center') 
    # Bond Dipole Moment (QM)
    Dipoles = PDB.get_bond_dipole(pdb_obj.qm_cluster_fchk, a1qm, a2qm)
    # write to csv or plot
    write_data(pdb_obj.MutaFlags, {'E': Es, 'Bond Dipole': Dipoles},  data_output_path)

    endtime_4 = datetime.datetime.now()
    print(f'Analysis: {endtime_4 - endtime_3}')

    endtime = datetime.datetime.now()
    print(f'total time: {endtime - starttime}')

if __name__ == "__main__":
    main()