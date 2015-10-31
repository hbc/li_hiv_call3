#!/bin/bash
set -eu -o pipefail

codedir=~/hsph/li_hiv_call3
python=/usr/local/share/bcbio_nextgen/anaconda/bin/python

mkdir -p variant_refs
cd variant_refs
ln -s -f ../li_hiv_call3
#$python $codedir/scripts/prepare_controls.py $codedir/inputs/control_libraries.csv

$python $codedir/scripts/control_to_vcf.py $codedir/config/reconstruct.yaml
