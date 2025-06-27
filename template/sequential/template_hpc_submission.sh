#!/bin/bash
#SBATCH --job-name=qkst_template
#SBATCH --account=yang_lab
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=2-00:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE

source ~/setup-accre.sh

python -u template_main.py > template_main.py.out 2>&1
