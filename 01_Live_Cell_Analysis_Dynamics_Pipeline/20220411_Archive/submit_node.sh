#!/bin/bash -l

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J 20211210
#SBATCH --mail-type=ALL
#SBATCH --mail-user=slurm@deliz.me
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --time=24:00:00

module purge
module load jdk/8.265 gcc/10 impi/2021.2 fftw-mpi R/4.0.2
echo 'modules loaded'

conda activate dynamics_pipeline
echo 'conda activated'

cd /raven/u/deliz/dynamics_pipeline

export OMP_NUM_THREDS=72

python mission_control.py '/raven/u/deliz/new_pipeline/20211210/Input/parameter_tables' 12