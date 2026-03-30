#!/bin/bash
#SBATCH --job-name=shr_main
#SBATCH --account=yang_lab
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --time=10-00:00:00
#SBATCH --signal=B:USR1@30
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --mail-user=qianzhen.shao@vanderbilt.edu
#SBATCH --mail-type=END,FAIL

# Gaussian (for formchk)
module load Gaussian/16.B.01
# AmberTool (for tleap antechamber parmchk)
source /home/shaoq1/bin/amber_env/amber-accre.sh
# Multiwfn
export KMP_STACKSIZE=200000000
export Multiwfnpath=/home/shaoq1/bin/Multiwfn_3.7_bin_Linux_noGUI
export PATH=$PATH:$Multiwfnpath
# EnzyHTP
source ~/bin/miniconda3/bin/activate new_enzy_htp

exec python -u enzyhtp_main.py > shrapnel_main.py.out 2>&1