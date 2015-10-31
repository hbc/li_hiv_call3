#!/bin/bash
set -eu -o pipefail

codedir=~/hsph/li_hiv_call3
tdir=therapyedge/Harv-3.22
python=/usr/local/share/bcbio_nextgen/anaconda/bin/python

$python $codedir/scripts/reconstruct_populations.py $codedir/config/reconstruct.yaml $tdir/2015-10-09_09-55_Control-Library-6480_NT.vcf $tdir/regions.bed
