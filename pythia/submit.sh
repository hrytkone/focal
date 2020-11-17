#!/bin/bash

if [ -z "$1" ]
then
    echo "Usage: ./`basename $0` comment njobs[=1] nevents[=1000]"
    exit 0
fi

if [ -z "$2" ]
then
    njobs=1
else
    njobs=$2
fi

if [ -z "$3" ]
then
    nevents=1000
else
    nevents=$3
fi


for (( i=1; i<=$njobs; i++ ))
do
    outputdir=run_${1}_job$i
    mkdir $outputdir
    mkdir ${outputdir}/logs
    sbatch -o ${outputdir}/logs/log -e ${outputdir}/logs/errout -J pionPion -n 1 run $i $nevents $outputdir
    sleep 1
done
