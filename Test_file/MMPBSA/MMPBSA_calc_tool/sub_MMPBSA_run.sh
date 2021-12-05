#!/bin/bash
#SBATCH --partition=production
#SBATCH --job-name=MMPBSA_abag
#SBATCH --nodes=1             
#SBATCH --tasks-per-node=24
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00 
#SBATCH --account=yang_lab_csb


############################################################################
#                                                                          #
######################### sbatch script for MD job #########################
#                                                                          #
############################################################################
# Write no bash commands above this hash box!

module load GCC/6.4.0-2.28  OpenMPI/2.1.1
module load Amber/17-Python-2.7.14
#module unload Python/2.7.14 numpy/1.13.1-Python-2.7.14

echo ACCRE Amber18 MMPBSA calculation of MMPBSA_abag on node: `hostname -s`
echo --------------------------------------------- 
echo The job started on: $(date) 
STARTTIME=$(date +%s)

python -u run_MMPBSA.py > run_MMPBSA.py.out

echo The job finished on: $(date)
ENDTIME=$(date +%s)
echo Overall the job took $(($ENDTIME - $STARTTIME)) seconds.
echo --------------------------------------------- 

################### Job Ended ###################
exit 0
