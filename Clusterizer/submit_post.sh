#!/bin/bash
## sbatch usage: sbatch submit.sh <input> <jobstart> <jobend>
## sbatch will check arguments from the comments in the
## beginning of this file.
#SBATCH --job-name=post
#SBATCH --account=project_2003583
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL #uncomment to enable mail
# If you change output/error here, please change
# the mv command at the end of this macro
#SBATCH --output=logs/output_run%a.txt
#SBATCH --error=logs/errors_run%a.txt

comment=$1  # simulation directory name
jobstart=$2 # first job to include
jobend=$3   # last job to include

simdir=/scratch/project_2003583/focal-full-sim/output/${comment}
clusterdir=/scratch/project_2003583/focal-full-sim/output_clustering/${comment}
outdir=/scratch/project_2003583/focal-full-sim/HeidiPostProcessing/output/${comment}

echo "Starting Slurm array task ${SLURM_ARRAY_TASK_ID}"

mkdir -p $outputdir
mkdir -p ${outputdir}/logs
/scratch/project_2003583/focal-full-sim/HeidiPostProcessing/run_post $simdir $clusterdir $outdir $jobstart $jobend
sleep 1
mv logs/output_run${SLURM_ARRAY_TASK_ID}.txt ${outputdir}/logs/.
mv logs/errors_run${SLURM_ARRAY_TASK_ID}.txt ${outputdir}/logs/.
sleep 1
