#!/bin/bash
#SBATCH --job-name=JM_FRA
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=pascal
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks=8
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00

export PGI_FASTMATH_CPU=sandybridge
export OMP_THREAD_LIMIT=256

# Gaussian
module load Gaussian/16.B.01
mkdir $TMPDIR/$SLURM_JOB_ID
export GAUSS_SCRDIR=$TMPDIR/$SLURM_JOB_ID
# Amber GPU
export AMBERHOME2=/dors/csb/apps/amber19/
export CUDA_HOME=$AMBERHOME2/cuda/10.0.130
export LD_LIBRARY_PATH=$AMBERHOME2/cuda/10.0.130/lib64:$AMBERHOME2/cuda/RHEL7/10.0.130/lib:$LD_LIBRARY_PATH
# Amber CPU
module load GCC/6.4.0-2.28  OpenMPI/2.1.1
module load Amber/17-Python-2.7.14
module unload Python/2.7.14 numpy/1.13.1-Python-2.7.14
# Multiwfn
export KMP_STACKSIZE=200000000
export Multiwfnpath=/home/jiany37/Software_installed/Multiwfn_3.7_bin_Linux_noGUI
export PATH=$PATH:$Multiwfnpath
# EnzyHTP
source ~/bin/miniconda3/bin/activate
conda activate MutaGen
export PYTHONPATH=$PYTHONPATH:~/bin/EnzyHTP

python -u FAcD-gen-fra.py > FAcD-gen-fra.py.out

rm -rf $TMPDIR/$SLURM_JOB_ID