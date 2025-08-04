#!/bin/bash -l

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J ImageAnalysis
#SBATCH --mail-type=none
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --time=24:00:00

module purge
module load jdk/8.265
module load gcc/10 impi/2021.2
conda activate r_env

cd /raven/u/deliz/dynamics_pipeline/r_scripts

Rscript --vanilla extract_intensity.R /raven/u/deliz/Input/parameter_tables

Rscript --vanilla colocalization_intensity-based.R /raven/u/deliz/Input/parameter_tables

Rscript --vanilla colocalization_coordinate-based.R /raven/u/deliz/Input/parameter_tables

Rscript --vanilla changepoint.R /raven/u/deliz/Input/parameter_tables

Rscript --vanilla compile_tables.R /raven/u/deliz/Input/parameter_tables
