#!/bin/bash
#SBATCH --partition=debug
#SBATCH --job-name=EnzyHTP-test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=30:00
#SBATCH --account=yang_lab_csb

source ~/bin/miniconda3/bin/activate
conda activate MutaGen
export PYTHONPATH=$PYTHONPATH:~/bin/EnzyHTP

pytest -s -m 'accre or clean' > ./test/pytest.out