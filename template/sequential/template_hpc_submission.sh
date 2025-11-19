#!/bin/bash
#SBATCH --job-name=qkst_template
#SBATCH --account=xxx
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=10-00:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE

# Gaussian (for formchk)
module load Gaussian/16.B.01
# AmberTool (for tleap antechamber parmchk)
source /home/shaoq1/bin/amber_env/amber-accre.sh
# Multiwfn
export KMP_STACKSIZE=200000000
export Multiwfnpath=/home/shaoq1/bin/Multiwfn_3.8_dev_bin_Linux_noGUI
export PATH=$PATH:$Multiwfnpath
# EnzyHTP
source ~/bin/miniconda3/bin/activate enzy_htp

python -u template_main.py > template_main.py.out 2>&1
