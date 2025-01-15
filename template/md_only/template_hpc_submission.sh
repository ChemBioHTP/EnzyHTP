#!/bin/bash
#SBATCH --job-name=qkst_template
#SBATCH --account=xxx # CHANGE this to your account
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
source /sb/apps/amber22/amber.sh
# Multiwfn
module load GCC/6.4.0-2.28  OpenMPI/2.1.1
export KMP_STACKSIZE=200000000
export Multiwfnpath=/home/shaoq1/bin/Multiwfn_3.7_bin_Linux_noGUI
export PATH=$PATH:$Multiwfnpath
# EnzyHTP
source ~/bin/miniconda3/bin/activate enzy_htp # CHANGE this to your command that activate the environment. (Note: if you are using conda from module load. You would need to copy the conda init section from ~/.bashrc to here before running conda activate xxx)

python -u template_main.py > template_main.py.out 2>&1
