import sys
import os
import glob
import datetime
import numpy as np
from Class_PDB import *
from Class_Conf import *
from Class_ONIOM_Frame import *
from helper import write_data, line_feed

# settings
Config.n_cores = 16
Config.max_core = 2000
Config.PC_cmd = 'srun'
Config.Amber.conf_equi['nstlim'] = 50000
Config.Amber.conf_prod['nstlim'] = 500000
Config.debug = 1
wkflow_log_path = './FAcD-Gen.py.log'

def main():
	starttime = datetime.datetime.now()

	of = open(wkflow_log_path, 'wb', buffering=0)
	for i in range(25):
		of.write(('Round: '+str(i)+line_feed).encode('utf-8'))												# Remove water and ion in original crystal file 											(self.path) 
		# --- Preparation ---
		pdb_obj = PDB('./FAcD-FA-ASP.pdb', wk_dir='./FAcD_test'+str(i))		# Initiate with the file from protein data bank 											(self.path) 
		pdb_obj.rm_wat()
		of.write(('Preparation: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))													# Remove water and ion in original crystal file 											(self.path) 
		# pdb_obj.get_protonation()											# For most crystal files, add hydrogens														(self.path) 

		# --- Sample with MD ---
		pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1) 										# Generate parameter files *savepdb save the exact structure use in MD for future analysis 	(self.path) 
		pdb_obj.PDBMD(tag=这个是文件夹的tag, engine='Amber_pmemd_gpu', equi_cpu=1)				# Run MD 																					(self.nc) 
		of.write(('MD: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))													# Remove water and ion in original crystal file 											(self.path) 
		# sample
		pdb_obj.nc2mdcrd(point=100)											# Sample from trajactory 																	(self.mdcrd) 
		
	of.close()




if __name__ == "__main__":
	main()
