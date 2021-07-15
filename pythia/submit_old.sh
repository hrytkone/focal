#!/bin/bash

if [ -z "$1" ]
then
    echo "Usage: ./`basename $0` comment njobs[=1] bUseLeading[=1]"
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
    buseleading=1
else
    buseleading=$3
fi

outputdir=output/pythiaOut_${1}
mkdir $outputdir

for (( i=1; i<=$njobs; i++ ))
do
    mkdir ${outputdir}/logs_job${i}
    sbatch --exclusive=user -o ${outputdir}/logs_job${i}/log -e ${outputdir}/logs_job${i}/errout -J pythiaf -n 1 run $i $outputdir $buseleading
    sleep 1
done
