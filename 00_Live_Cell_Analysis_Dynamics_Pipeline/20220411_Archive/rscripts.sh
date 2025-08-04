#!/bin/bash -l

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J ImgPipe
#SBATCH --mail-type=END
#SBATCH --mail-user=deliz@mpiib-berlin.mpg.de
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=144
#SBATCH --time=24:00:00
#SBATCH --mem=500000

module purge
module load jdk/8.265 gcc/10 impi/2021.2 fftw-mpi R/4.0.2
echo 'modules loaded'

conda activate dynamics_pipeline
echo 'conda activated'

cd /raven/u/deliz/dynamics_pipeline

export OMP_NUM_THREDS=144

Rscript --vanilla --verbose /raven/u/deliz/dynamics_pipeline/r_scripts/colocalization_intensity-based.R /raven/u/deliz/new_pipeline/Input/parameter_tables
Rscript --vanilla --verbose /raven/u/deliz/dynamics_pipeline/r_scripts/colocalization_coordinate-based.R /raven/u/deliz/new_pipeline/Input/parameter_tables
Rscript --vanilla --verbose /raven/u/deliz/dynamics_pipeline/r_scripts/changepoint.R /raven/u/deliz/new_pipeline/Input/parameter_tables
Rscript --vanilla --verbose /raven/u/deliz/dynamics_pipeline/r_scripts/compile_tables.R /raven/u/deliz/new_pipeline/Input/parameter_tables

