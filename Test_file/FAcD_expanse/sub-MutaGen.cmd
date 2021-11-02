#!/bin/bash
#SBATCH --job-name=FAcD_gpu_Test          # Assign an 8-character name to your job
#SBATCH --account=cla296
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gpus=1
#SBATCH --mem=50G
#SBATCH --no-requeue
#SBATCH --time=18:00:00         # Total run time limit (HH:MM:SS)


module purge
module load gpu  
module load openmpi
module load slurm
module load amber/20
module load gaussian/16.C.01-cuda
export TMPDIR=/scratch/$USER/job_$SLURM_JOB_ID
mkdir $TMPDIR
export GAUSS_SCRDIR=$TMPDIR

source ~/bin/miniconda3/bin/activate
conda activate MutaGen

export PYTHONPATH=$PYTHONPATH:~/enzyme_workflow

python -u FAcD-Gen-fixlog.py >> FAcD-Gen-fixlog.py.out
