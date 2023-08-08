from core.clusters.accre import Accre
import datetime
from Class_PDB import PDB
from Class_Conf import Config
from Class_ONIOM_Frame import *
from helper import write_data

# settings
Config.n_cores = 8
Config.max_core = 2000
Config.Amber.AmberEXE_GPU = '/dors/csb/apps/amber19/bin/pmemd.cuda'
Config.Amber.conf_equi['nstlim'] = 500000
Config.Amber.conf_prod['nstlim'] = 50000000
Config.debug = 1
data_output_path = './Mutation-E-BD.dat'
time_log_path = './FAcD_test_timelog.log'

def main():
    starttime = datetime.datetime.now() ##
    of = open(time_log_path, 'wb', buffering=0) ##
    # --- Preparation ---
    pdb_obj = PDB('./FAcD-FA-ASP.pdb')
    of.write(('Preparation: PDB: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8')) ##
    lasttime = datetime.datetime.now() ## 

    # --- Operation ---
    # Mutation
    Muta_tag = pdb_obj.Add_MutaFlag('K83D')
    pdb_obj.PDB2PDBwLeap()
    time_2 = datetime.datetime.now()
    of.write(('Mutation: p2pwl: '+str(datetime.datetime.now()- lasttime)+line_feed).encode('utf-8')) ##
    lasttime = datetime.datetime.now() ##
    # protonation modification
    pdb_obj.rm_allH()
    pdb_obj.get_protonation(if_prt_ligand=0)
    # use minimization to relax each mutated PDB
    pdb_obj.PDB2FF(local_lig=0)
    pdb_obj.PDBMin(engine='Amber_GPU', if_cluster_job=0)
    pdb_obj.rm_wat()
    of.write(('Mutation: PDBMin: '+str(datetime.datetime.now()- lasttime)+line_feed).encode('utf-8')) ##
    lasttime = datetime.datetime.now() ##

    # --- Sample with MD ---
    pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1)
    pdb_obj.PDBMD(engine='Amber_GPU', equi_cpu=1, if_cluster_job=0)
    of.write(('MD: '+str(datetime.datetime.now()- lasttime)+line_feed).encode('utf-8')) ##
    lasttime = datetime.datetime.now() ##
    
    # --- QM cluster ---
    pdb_obj.nc2mdcrd(point=100)
    atom_mask = ':108,298'
    g_route = '# pbe1pbe/def2TZVP nosymm'
    pdb_obj.PDB2QMCluster(atom_mask, g_route=g_route, ifchk=1, if_cluster_job=0)
    pdb_obj.get_fchk(keep_chk=0)
    of.write(('QM Cluster: '+str(datetime.datetime.now()- lasttime)+line_feed).encode('utf-8')) ##
    lasttime = datetime.datetime.now() ##

    # --- Analysis ---
    # targeting C-F bond
    a1 = int(pdb_obj.stru.ligands[0].CH3)
    a2 = int(pdb_obj.stru.ligands[0].F)
    a1qm = pdb_obj.qm_cluster_map[str(a1)]
    a2qm = pdb_obj.qm_cluster_map[str(a2)]
    # Field Strength (MM)
    E_atom_mask = ':1-107,109-297'
    Es = pdb_obj.get_field_strength(E_atom_mask, a1=a1 ,a2=a2 ,bond_p1='center')
    of.write(('Analysis: Get field strength: '+str(datetime.datetime.now()- lasttime)+line_feed).encode('utf-8')) ##
    lasttime = datetime.datetime.now() ##
    # Bond Dipole Moment (QM)
    Dipoles = PDB.get_bond_dipole(pdb_obj.qm_cluster_fchk, a1qm, a2qm)	# Run Bond Dipole Moment analysis
    of.write(('Get bond dipole: '+str(datetime.datetime.now()- lasttime)+line_feed).encode('utf-8')) ##
    lasttime = datetime.datetime.now() ##
    # Mutation distance
    r1 = pdb_obj.stru.ligands[0]
    r2 = pdb_obj.stru.chains[ord(pdb_obj.MutaFlags[0][1])-65][int(pdb_obj.MutaFlags[0][2])-1]
    Dist = pdb_obj.stru.get_resi_dist(r1, r2)
    # write to csv or plot
    write_data(pdb_obj.MutaFlags, {'E': Es, 'Bond Dipole': Dipoles, 'Distance': Dist}, data_output_path)
    		
    endtime = datetime.datetime.now() ##
    print(endtime - starttime) ##
    of.close() ##

if __name__ == "__main__":
    main()
