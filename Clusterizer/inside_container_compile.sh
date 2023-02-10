#!/bin/bash

eval `alienv -w /home/hadi/alicesw/sw --no-refresh printenv FOCAL/latest,AliDPG/latest`

set -x
cmd="root -l -b -q 'CompileEvent.C()'"
eval $cmd