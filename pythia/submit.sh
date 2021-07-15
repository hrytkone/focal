#!/bin/bash
## usage : sbatch submit.sh <comment> <buseleading>
#SBATCH --job-name=focal
#SBATCH --output=output/logs/output_%j.txt
#SBATCH --error=output/logs/error_%j.txt
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000
#SBATCH --array=1-10
#SBATCH --time=10:00:00

if [ -z "$1" ]
then
    echo "Usage: ./`basename $0` comment bUseLeading[=1]"
    exit 0
fi

if [ -z "$2" ]
then
    buseleading=1
else
    buseleading=$2
fi

outputdir=output/pythiaOut_${1}
mkdir -p $outputdir

n=$SLURM_ARRAY_TASK_ID
/n/work00/heimarry/focal/pythia/run $n $outputdir $buseleading
sleep 1
