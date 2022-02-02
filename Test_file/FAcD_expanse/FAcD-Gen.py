import sys
import os
import glob
import datetime
from Class_PDB import *
from Class_Conf import *
from Class_ONIOM_Frame import *
from helper import write_data

# settings
Config.n_cores = 16
Config.max_core = 2000
Config.PC_cmd = 'srun'
Config.Amber.conf_equi['nstlim'] = 50000
Config.Amber.conf_prod['nstlim'] = 500000
Config.debug = 1
data_output_path = './Mutation-E-BD.dat'

def main():
	starttime = datetime.datetime.now()

	for i in range(1):
		# --- Preparation ---
		pdb_obj = PDB('./FAcD-FA-ASP.pdb', wk_dir='./FAcD_test'+str(i))		# Initiate with the file from protein data bank 											(self.path) 
		pdb_obj.rm_wat()													# Remove water and ion in original crystal file 											(self.path) 
		# pdb_obj.get_protonation()											# For most crystal files, add hydrogens														(self.path) 

		# --- Operation ---
		# Mutation
		Muta_tag = pdb_obj.Add_MutaFlag('r')								# Generate a target "Flag" for mutation 													
		if pdb_obj.MutaFlags[0][2] == '108': 								# Keep the key residue 
			continue
		pdb_obj.PDB2PDBwLeap()												# Deploy mutation 																			(self.path) 
		# protonation modification
		pdb_obj.rm_allH()													# Remove all hydrogens after mutation (residues only) 										(self.path) 
		pdb_obj.get_protonation()											# Determine protonation state again 														(self.path) 
		# use minimization to relax each mutated PDB
		pdb_obj.PDB2FF()													# Generate parameter files for MD simulation 
		pdb_obj.PDBMin(engine='Amber_pmemd_gpu')							# Minimization 																				(self.path) 

		# --- Sample with MD ---
		pdb_obj.rm_wat()													# Remove water from the minimization 														(self.path) 												
		pdb_obj.PDB2FF(ifsavepdb=1) 										# Generate parameter files *savepdb save the exact structure use in MD for future analysis 	(self.path) 
		pdb_obj.PDBMD(tag=Muta_tag, engine='Amber_pmemd_gpu', equi_cpu=1)				# Run MD 																					(self.nc) 
		# sample
		pdb_obj.nc2mdcrd(point=100)											# Sample from trajactory 																	(self.mdcrd) 

		# --- QM cluster ---
		atom_mask = ':108,298'												# Define QM cluster / can also use some presets: ligand; residues within a distance using a Layer object
		g_route = '# hf/6-31G(d) pop=cm5'									# QM keywords
		pdb_obj.PDB2QMCluster(atom_mask, g_route=g_route, ifchk=1)			# Run QM cluster calculation																(self.qm_cluster_out, self.qm_cluster_chk) 
		pdb_obj.get_fchk(keep_chk=0)										# Save fchk files for analysis																(self.qm_cluster_fchk) 
		
		# --- Analysis ---
		# targeting C-F bond
		a1 = int(pdb_obj.stru.ligands[0].CH3)
		a2 = int(pdb_obj.stru.ligands[0].F)
		a1qm = pdb_obj.qm_cluster_map[str(a1)]
		a2qm = pdb_obj.qm_cluster_map[str(a2)]
		# Field Strength (MM)
		E_atom_mask = ':1-107,109-297'										# Define atoms for field strength calculation
		Es = pdb_obj.get_field_strength(E_atom_mask, a1=a1 ,a2=a2 ,bond_p1='center') # Run Field Strength analysis
		# Bond Dipole Moment (QM)
		Dipoles = PDB.get_bond_dipole(pdb_obj.qm_cluster_fchk, a1qm, a2qm)	# Run Bond Dipole Moment analysis
		# Mutation distance
		r1 = pdb_obj.stru.ligands[0]
		r2 = pdb_obj.stru.chains[ord(pdb_obj.MutaFlags[0][1])-65][int(pdb_obj.MutaFlags[0][2])-1]
		Dist = pdb_obj.stru.get_resi_dist(r1, r2)

		# write to csv or plot												
		write_data(pdb_obj.MutaFlags, {'E': Es, 'Bond Dipole': Dipoles, 'Distance': Dist}, data_output_path)		# Current data: Mutation - MD geometry - QM cluster wavefunction = Field strength at bond - Bond dipole moment
	
	endtime = datetime.datetime.now()
	print(endtime - starttime)



if __name__ == "__main__":
	main()
