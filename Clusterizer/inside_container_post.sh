#!/bin/bash

eval `alienv -w /home/hadi/alicesw/sw --no-refresh printenv FOCAL/latest,AliDPG/latest`

set -x
cmd="aliroot -l -b -q 'RunCombination.C(\"${SIMDIR}\",\"${CLUSTERDIR}\",\"${OUTDIR}\",'${JOBSTART}','${JOBEND}')'"
eval $cmd

