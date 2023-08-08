import sys
import os
import glob
import datetime
from Class_PDB import *
from Class_Conf import *
from Class_ONIOM_Frame import *

Config.n_cores = 24
Config.max_core = 2000
#Config.PC_cmd = 'srun' # srun does not work like this use mpi instead
#Config.Amber.AmberHome = '$AMBERHOME2'
#Config.Amber.AmberEXE = '/dors/csb/apps/amber19/bin/pmemd.cuda_DPFP.MPI'
Config.Amber.conf_prod['nstlim'] = 500000
Config.debug = 1


def main():
	starttime = datetime.datetime.now()

	for i in range(1):
		pdb_obj = PDB('./FAcD-FA-ASP_rmW_HA227C_aH_ff_min_rmW.pdb')
		# pdb_obj.rm_wat()
		#Mutation
		#Muta_tag = pdb_obj.Add_MutaFlag('r')
		#print(Muta_tag)
		#if pdb_obj.MutaFlags[0][2] == '108':
		#	continue
		#pdb_obj.PDB2PDBwLeap()
		# protonation modification
		#pdb_obj.rm_allH()
		#pdb_obj.get_protonation()
		#use minimization to relax each mutated PDB
		#pdb_obj.PDB2FF(ifsavepdb=1)
		#pdb_obj.PDBMin(engine='Amber_pmemd_cpu')

		# #run MD
		#pdb_obj.rm_wat()
		pdb_obj.PDB2FF(ifsavepdb=1)
		#pdb_obj.PDBMD(tag=Muta_tag, engine='Amber_pmemd_cpu')

		# sample
		# pdb_obj.mdcrd='../FAcD_test0/MD_HA227C/prod.mdcrd'

		# QMcluster
		pdb_obj.get_stru()
		atom_mask = ':108,298'
		sele_list, sele_map = pdb_obj.stru.get_sele_list(atom_mask, fix_end='H', prepi_path=pdb_obj.prepi_path)
		pdb_obj.qm_cluster_map = sele_map	
		# g_route = '# hf/6-31G(d) pop=cm5'
		# pdb_obj.PDB2QMCluster(atom_mask, g_route=g_route, ifchk=1)

		# get value
		# pdb_obj.qm_cluster_fchk = glob.glob('./QM_cluster/qm_cluster_0.fchk')
		# values = pdb_obj.get_field_strength()
		print(len(pdb_obj.stru.ligands[0]))
		a1 = int(pdb_obj.stru.ligands[0].CH3)
		a2 = int(pdb_obj.stru.ligands[0].F)
		a1qm = pdb_obj.qm_cluster_map[str(a1)]
		a2qm = pdb_obj.qm_cluster_map[str(a2)]		
		values2 = pdb_obj.get_bond_dipole(['./QM_cluster/qm_cluster_0.fchk'], a1qm, a2qm)
		print(values2)		
		# write to csv or plot
		# xxx
	
	endtime = datetime.datetime.now()
	print(endtime - starttime)



if __name__ == "__main__":
	main()
