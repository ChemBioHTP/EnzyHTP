#!/bin/bash
#SBATCH --partition=production
#SBATCH --job-name=MMPBSA_prep
#SBATCH --nodes=1             
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=10:00:00 
#SBATCH --account=yang_lab_csb


############################################################################
#                                                                          #
######################### sbatch script for MD job #########################
#                                                                          #
############################################################################
# Write no bash commands above this hash box!

module load GCC/6.4.0-2.28  OpenMPI/2.1.1
module load Amber/17-Python-2.7.14
module unload Python/2.7.14 numpy/1.13.1-Python-2.7.14

# EnzyHTP
source ~/bin/miniconda3/bin/activate
conda activate MutaGen
export PYTHONPATH=$PYTHONPATH:~/enzyme_workflow

echo ACCRE Amber20 MMPBSA calculation of MMPBSA_abag on node: `hostname -s`
echo --------------------------------------------- 
echo The job started on: $(date) 
STARTTIME=$(date +%s)

python -u MMPBSA_prep.py > MMPBSA_prep.py.out

echo The job finished on: $(date)
ENDTIME=$(date +%s)
echo Overall the job took $(($ENDTIME - $STARTTIME)) seconds.
echo --------------------------------------------- 

################### Job Ended ###################
exit 0
