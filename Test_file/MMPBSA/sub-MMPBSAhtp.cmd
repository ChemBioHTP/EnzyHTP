#!/bin/bash
#SBATCH --job-name=MMPBSA-htp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=pascal
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --mem=50G
#SBATCH -t 3-00:00:00
#SBATCH --export=ALL

export AMBERHOME2=/dors/csb/apps/amber19/
export CUDA_HOME=$AMBERHOME2/cuda/10.0.130
export LD_LIBRARY_PATH=$AMBERHOME2/cuda/10.0.130/lib64:$AMBERHOME2/cuda/RHEL7/10.0.130/lib:$LD_LIBRARY_PATH

module load Intel/2017.4.196  IntelMPI/2017.3.196
module load Amber/17-Python-2.7.14
module unload Python/2.7.14 numpy/1.13.1-Python-2.7.14

#EnzyHTP
source ~/bin/miniconda3/bin/activate
conda activate MutaGen
export PYTHONPATH=$PYTHONPATH:~/enzyme_workflow

python -u MMPBSA-MD-gen.py > MMPBSA-MD-gen.py.out
