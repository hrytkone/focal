#!/bin/bash
## sbatch usage: sbatch submit.sh <input> <jobstart> <jobend>
## sbatch will check arguments from the comments in the
## beginning of this file.
#SBATCH --job-name=compile
#SBATCH --account=project_2003583
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=FAIL #uncomment to enable mail
# If you change output/error here, please change
# the mv command at the end of this macro
#SBATCH --output=logs/output_compilation_%a.txt
#SBATCH --error=logs/errors_compilation_%a.txt

echo "Starting Slurm array task ${SLURM_ARRAY_TASK_ID}"

/scratch/project_2003583/focal-full-sim/HeidiPostProcessing/run_compilation
