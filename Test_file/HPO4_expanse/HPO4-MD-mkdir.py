from helper import mkdir
import os

def main():
	# get mutant list
	MutaFlags = []
	with open('Single_Mut_list.txt') as f:
		for line in f:
			MutaFlag = line.strip()
			if MutaFlag != '':
				MutaFlags.append(MutaFlag)
	# deploy mutation
	for Flag in MutaFlags:
		mkdir(Flag)
		with open('./'+Flag+'/HPO4-MD-gen.py','w') as of:
			of.write('''import datetime
from Class_PDB import *
from Class_Conf import *
from helper import line_feed

# settings
Config.n_cores = 1
Config.max_core = 2000
#Config.PC_cmd = 'srun' # srun does not work like this use mpi instead
wkflow_log_path = './MD-gen.log'

# MD settings
Config.Amber.box_type = 'oct'
Config.Amber.conf_min['ntr'] = '0'
Config.Amber.conf_heat['restraintmask'] = "':MOL'"
Config.Amber.conf_heat['restraint_wt'] = "100.0"
Config.Amber.conf_heat['tempi'] = '100.0'
Config.Amber.conf_heat['temp0'] = '296.15'

Config.Amber.conf_equi['ntr'] = '0'
Config.Amber.conf_equi['nstlim'] = 500000  # 1ns
Config.Amber.conf_equi['temp0'] = '296.15'

Config.Amber.conf_prod['nstlim'] = 30000000 # 60ns
Config.Amber.conf_prod['ntwx'] = '25000'
Config.Amber.conf_prod['temp0'] = '296.15'

def main():
        starttime = datetime.datetime.now()
        of = open(wkflow_log_path, 'wb', buffering=0)
        # deploy mutation
        of.write(('Working on: '+ '''+repr(Flag)+''' +line_feed).encode('utf-8'))
        # --- Preparation ---
        pdb_obj = PDB('WT-HPO4.pdb')
        of.write(('Preparation: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))

        # --- Operation ---
        # Mutation
        pdb_obj.Add_MutaFlag('''+repr(Flag)+''')
        pdb_obj.PDB2PDBwLeap()
        of.write(('Mutation: p2pwl: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))

        # use minimization to relax each mutated PDB
        pdb_obj.PDB2FF(local_lig=0)
        pdb_obj.PDBMin(engine='Amber_GPU')
        of.write(('Mutation: PDBMin: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))

        # --- Sample with MD ---
        pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1)
        pdb_obj.PDBMD(engine='Amber_GPU')
        of.write(('MD: '+str(datetime.datetime.now()- starttime)+line_feed).encode('utf-8'))			
			
        endtime = datetime.datetime.now()
        print(endtime - starttime)
        of.close()

if __name__ == "__main__":
        main()''')

		with open('./'+Flag+'/sub-EnzyHTP.cmd','w') as of:
			of.write('''#!/bin/bash
#SBATCH --job-name=MD_'''+Flag+'''          # Assign an 8-character name to your job
#SBATCH --account=cla296
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --mem=50G
#SBATCH --no-requeue
#SBATCH --time=24:00:00         # Total run time limit (HH:MM:SS)

module purge
module load gpu
module load openmpi
module load slurm
module load amber/20

#EnzyHTP
source ~/bin/miniconda3/bin/activate
conda activate MutaGen
export PYTHONPATH=$PYTHONPATH:~/enzyme_workflow

python -u HPO4-MD-gen.py > HPO4-MD-gen.py.out''')
		os.system('cp ./*pdb '+Flag)



if __name__ == "__main__":
	main()
