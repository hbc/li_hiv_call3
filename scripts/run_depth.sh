#!/bin/bash
set -eu -o pipefail

parallel -t -j 1 "bammarkduplicates I={} O={.}-nodup.bam rmdup=1" ::: Control*-ready.bam
parallel -t -j 1 "sambamba index {}" ::: Control*-ready-nodup.bam
parallel -t -j 1 "echo {}; sambamba depth window -w 200 --overlap 100 {}" ::: Control*-ready.bam Control*-ready-nodup.bam
