from glob import glob
import datetime
from Class_PDB import *
from Class_Conf import *
from Class_ONIOM_Frame import *
from helper import line_feed

# settings
Config.n_cores = 16
Config.max_core = 2000
Config.PC_cmd = 'srun'
wkflow_log_path = './MD-gen.log'
Config.Amber.AmberHome='$AMBERHOME2'

# MD settings
Config.Amber.box_type = 'oct'
Config.Amber.conf_min['ntr'] = '0'
Config.Amber.conf_heat['ntr'] = '0'
Config.Amber.conf_equi['ntr'] = '0'
Config.Amber.conf_equi['nstlim'] = 50000  # 0.1ns
Config.Amber.conf_prod['nstlim'] = 500000 # 1ns


def main():
	starttime = datetime.datetime.now()

	of = open(wkflow_log_path, 'wb', buffering=0)
	for inp_pdb in glob('*.pdb'):
		of.write(('Working on: '+inp_pdb+line_feed).encode('utf-8'))	
		# --- Preparation ---
		pdb_obj = PDB(inp_pdb, wk_dir='./'+inp_pdb[:-4])		
		pdb_obj.rm_wat()
		of.write(('Preparation: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))
		pdb_obj.rm_allH()
		pdb_obj.get_protonation()											

		# --- Sample with MD ---
		pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1)								
		pdb_obj.PDBMD(engine='Amber_pmemd_gpu')
		of.write(('MD: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))
		# sample
		pdb_obj.nc2mdcrd(point=100)
			
	endtime = datetime.datetime.now()
	print(endtime - starttime)
	of.close()




if __name__ == "__main__":
	main()
