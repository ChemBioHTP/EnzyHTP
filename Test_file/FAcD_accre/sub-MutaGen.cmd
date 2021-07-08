#!/bin/bash
#SBATCH --partition=production    # Partition (job queue)
#SBATCH --job-name=FAcDTest          # Assign an 8-character name to your job
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --tasks-per-node=24       # Number of tasks (usually = cores) on each node
#SBATCH --mem-per-cpu=4G   # Minimum memory (in MB) required per allocated CPU
#SBATCH --time=24:00:00         # Total run time limit (HH:MM:SS)
#SBATCH --export=ALL              # Export you current env to the job env
#SBATCH --account=yang_lab_csb

module load GCC/6.4.0-2.28  OpenMPI/2.1.1
#module load Intel/2017.4.196  IntelMPI/2017.3.196
module load Amber/17-Python-2.7.14
module unload Python/2.7.14 numpy/1.13.1-Python-2.7.14
#export AMBERHOME2=/dors/csb/apps/amber19/
#export CUDA_HOME=$AMBERHOME2/cuda/10.0.130
#export LD_LIBRARY_PATH=$AMBERHOME2/cuda/10.0.130/lib64:$AMBERHOME2/cuda/RHEL7/10.0.130/lib:$LD_LIBRARY_PATH

module load Gaussian/16.B.01
export PGI_FASTMATH_CPU=sandybridge
export OMP_THREAD_LIMIT=256

source ~/bin/miniconda3/bin/activate
conda activate MutaGen

export PYTHONPATH=$PYTHONPATH:~/enzyme_workflow

python FAcD-Gen.py > FAcD-Gen.py.out
