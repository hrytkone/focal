#!/bin/bash

echo "Compile event class"

# Settings for singularity:
export APPTAINER_BIND="/projappl/project_2003583,/scratch/project_2003583,$TMPDIR"
export APPTAINER_SHELL="/bin/bash --norc"
export APPTAINER_CACHEDIR=$TMPDIR

apptainer exec --home $HOME --workdir $TMPDIR -B /scratch/project_2003583/focal-full-sim:/home/hadi -B /appl/:/appl/  /scratch/project_2003583/focal-full-sim/container_cs8.sif /scratch/project_2003583/focal-full-sim/HeidiPostProcessing/inside_container_compile.sh

echo "Done"
