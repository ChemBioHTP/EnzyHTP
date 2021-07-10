#!/bin/bash
#SBATCH --job-name=FAcD_gpu_Test          # Assign an 8-character name to your job
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=pascal
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --mem=25G
#SBATCH --time=24:00:00         # Total run time limit (HH:MM:SS)
#SBATCH --export=ALL              # Export you current env to the job env

module load GCC/6.4.0-2.28  OpenMPI/2.1.1
module load Amber/17-Python-2.7.14
module unload Python/2.7.14 numpy/1.13.1-Python-2.7.14
export AMBERHOME2=/dors/csb/apps/amber19/
export CUDA_HOME=$AMBERHOME2/cuda/10.0.130
export LD_LIBRARY_PATH=$AMBERHOME2/cuda/10.0.130/lib64:$AMBERHOME2/cuda/RHEL7/10.0.130/lib:$LD_LIBRARY_PATH

module load Gaussian/16.B.01
export PGI_FASTMATH_CPU=sandybridge
export OMP_THREAD_LIMIT=256

source ~/bin/miniconda3/bin/activate
conda activate MutaGen

export PYTHONPATH=$PYTHONPATH:~/enzyme_workflow

python FAcD-Gen.py > FAcD-Gen.py.out
