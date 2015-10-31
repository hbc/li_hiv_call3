#!/bin/bash
set -eu -o pipefail

codedir=~/hsph/li_hiv_call3
python=/usr/local/share/bcbio_nextgen/anaconda/bin/python

$python $codedir/scripts/reconstruct_populations.py 1 $codedir/config/reconstruct.yaml Control2-ready.bam Control1-ready.bam

# -- abra testing
cores=4
JAR=abra-0.94-SNAPSHOT-jar-with-dependencies.jar
ref=$codedir/inputs/hiv_hxb2.fa
regions=$codedir/inputs/regions.bed

#window_regions=regions/all-regions-windows.bed
#rm -rf abra_temp_dir
#bedtools makewindows -b $regions -w 200 > $window_regions
#test_bam=Control1-ready-abra.bam
#test_bam_sort=Control1-ready-abra-sort.bam
#[[ -f $test_bam ]] || java -Xmx4G -jar $JAR --in Control1-ready.bam --out $test_bam --ref $ref --targets $window_regions --threads $cores --working abra_temp_dir --mad 250 --adc 1000000 -mnf 50 --kmer 35 > abra.log 2>&1
#[[ -f $test_bam_sort ]] || sambamba sort -t $cores -o $test_bam_sort $test_bam
#[[ -f $test_bam_sort.bai ]] || sambamba index -t $cores $test_bam_sort
#$python $codedir/scripts/reconstruct_populations.py 1 $codedir/config/reconstruct.yaml $test_bam_sort > validate-abra.txt
