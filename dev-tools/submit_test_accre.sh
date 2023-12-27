#!/bin/bash
#SBATCH --partition=debug
#SBATCH --job-name=devtest_enzy_htp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=30:00
#SBATCH --account=yang_lab_csb

# ! NOTE ! dont change this file. copy this file and use it. This way it wont mess up the repo

source /home/shaoq1/bin/amber_env/amber22.sh
source ~/bin/miniconda3/bin/activate new_EnzyHTP

python -u -m pytest test/_interface/test_amber_interface.py::test_amber_md_step_run -s > ./test/pytest.out
