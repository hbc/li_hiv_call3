#!/usr/bin/env python
"""Reconstruct HIV populations from aligned BAM files in regions of interest.

Requires:
  bedtools
  samtools
  bfc (https://github.com/lh3/bfc)
  pear (http://www.exelixis-lab.org/web/software/pear)
  bwa
  ViQuaS (unpack in work directory)
         (wget 'http://downloads.sourceforge.net/project/viquas/ViQuaS1.3.tar.gz' && tar -xzvpf)
  R
  perl

Library installs:
  cpanm -i BioPerl
  cpanm -i Perl4::CoreLibs # getopts.pl dependency
  R -e 'install.packages("seqinr")'
  R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Biostrings")'

  # local installs
  export PERL5LIB=/path/to/lib/perl5:
  export R_LIBS=/path/to/lib/R/site-library

  cpanm -i --local-lib=/path/to ...
  R -e '.libPaths(c("/path/to/lib/R/site-library")); ...'

Usage:
    reconstruct_populations.py <num_cores> <config_file> [<BAM files>]
"""
import contextlib
import os
import shutil
import sys
import subprocess

import joblib
import pyfaidx
import yaml

def main(cores, config_file, *bam_files):
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    fai_file = _index_ref_file(config["ref_file"])
    regions = _subset_region_file(config["regions"])
    to_run = []
    for bam_file in bam_files:
        for region, region_file in regions:
            subbam_file = _select_regions(bam_file, region_file, region)
            fq1, fq2 = _prep_fastq(subbam_file)
            if config["params"]["prep"] == "pear":
                prep_fq1 = _merge_fastq(fq1, fq2, subbam_file)
                prep_fq2 = None
            elif config["params"]["prep"] == "bfc":
                prep_fq1 = _correct_fastq(fq1, fai_file)
                prep_fq2 = _correct_fastq(fq2, fai_file)
            merged_bam_file = _realign_merged(prep_fq1, prep_fq2, config["ref_file"])
            to_run.append((merged_bam_file, config, region))
    def by_size((f, c, r)):
        chrom, start, end = r.split("-")
        return int(end) - int(start)
    to_run.sort(key=by_size)
    for recon_file in joblib.Parallel(int(cores))(joblib.delayed(_run_viquas)(f, c, r) for f, c, r in to_run):
        print recon_file

def _run_viquas(bam_file, config, region):
    out_dir = os.path.join(os.getcwd(), "%s-viquas" % os.path.splitext(bam_file)[0])
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file = os.path.join(out_dir, "ViQuaS-Spectrum.fa")
    bam_file = os.path.join(os.getcwd(), bam_file)
    viquas_dir = os.path.join(os.getcwd(), config["viquas_dir"])
    ref_file = os.path.join(os.getcwd(), config["ref_file"])
    # parameter settings for Illumina data from ViQuaS paper
    o, r = 5, 0.7
    # more lenient parameters to improve run time
    o, r = 3, 0.8
    chrom, start, end = region.split("-")
    size = int(end) - int(start)
    if not os.path.exists(out_file):
        with _chdir(out_dir):
            with _copy_viquas(viquas_dir):
                cmd = "Rscript ViQuaS.R {ref_file} {bam_file} {o} '' 0 {size}"
                subprocess.check_call(cmd.format(**locals()), shell=True)
    return _subset_viquas_file(out_file, region, config["ref_file"])

@contextlib.contextmanager
def _chdir(new_dir):
    cur_dir = os.getcwd()
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)

@contextlib.contextmanager
def _copy_viquas(viquas_dir):
    """ViQuaS requires running within the install directory.
    This copies it into the current directory so multiple instances can be
    run simultaneously.
    """
    files = ["ViQuaS.R"]
    dirs = ["tools", "viquas_files"]
    for fname in files:
        if not os.path.exists(fname):
            shutil.copy(os.path.join(viquas_dir, fname), os.path.join(os.getcwd(), fname))
    for dirname in dirs:
        if not os.path.exists(dirname):
            shutil.copytree(os.path.join(viquas_dir, dirname),
                            os.path.join(os.getcwd(), dirname))
    try:
        yield
    finally:
        for fname in files:
            if os.path.exists(fname):
                os.remove(fname)
        for dirname in dirs:
            if os.path.exists(dirname):
                shutil.rmtree(dirname)

def _prep_fastq(bam_file):
    """Extract paired end reads.
    """
    fq1 = "%s-1.fastq" % os.path.splitext(bam_file)[0]
    fq2 = "%s-2.fastq" % os.path.splitext(bam_file)[0]
    if not os.path.exists(fq1) or not os.path.exists(fq2):
        cmd = ("samtools sort -n -O bam -T {fq1} {bam_file} | "
               "bedtools bamtofastq -i /dev/stdin -fq {fq1} -fq2 {fq2}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return fq1, fq2

def _correct_fastq(in_file, fai_file):
    """Error correct input fastq file with bfc.
    """
    gsize = 0
    with open(fai_file) as in_handle:
        for line in in_handle:
            gsize += int(line.split()[1])
    out_file = "%s-corrected%s" % os.path.splitext(in_file)
    if not os.path.exists(out_file):
        cmd = "bfc -s {gsize} {in_file} > {out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def _merge_fastq(fq1, fq2, bam_file):
    """Extract paired end reads and merge using pear.
    """
    merged_out = "%s.assembled.fastq" % os.path.splitext(bam_file)[0]
    if not os.path.exists(merged_out):
        cmd = ["pear", "-f", fq1, "-r", fq2, "-o", os.path.splitext(bam_file)[0]]
        subprocess.check_call(cmd)
    return merged_out

def _realign_merged(fq1_file, fq2_file, ref_file):
    """Realign merged reads back to reference genome.
    """
    out_file = "%s.bam" % os.path.splitext(fq1_file)[0]
    if not fq2_file:
        fq2_file = ""
    if not os.path.exists(out_file):
        cmd = "bwa mem {ref_file} {fq1_file} {fq2_file} | samtools sort - -o {out_file} -O bam -T {out_file}-tmp"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def _select_regions(bam_file, region_file, region):
    """Extracts reads only in regions of interest using bedtools intersect.
    """
    base, ext = os.path.splitext(bam_file)
    out_file = "%s-%s%s" % (base, region, ext)
    if not os.path.exists(out_file):
        cmd = "bedtools intersect -abam {bam_file} -b {region_file} -f 1.0 > {out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def _subset_region_file(in_file):
    """Subset a BED file into individual regions to be analyzed separately.
    """
    out_dir = os.path.join(os.getcwd(), "regions")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out = []
    with open(in_file) as in_handle:
        for line in in_handle:
            chrom, start, end = line.strip().split("\t")[:3]
            base, ext = os.path.splitext(os.path.basename(in_file))
            out_file = os.path.join(out_dir, "%s-%s-%s-%s%s" % (base, chrom, start, end, ext))
            with open(out_file, "w") as out_handle:
                out_handle.write(line)
            out.append(("%s-%s-%s" % (chrom, start, end), out_file))
    return out

# reference management

def _subset_viquas_file(in_file, region, ref_file):
    """Subset output file by regions that differ from the reference genome.
    """
    base, ext = os.path.splitext(in_file)
    out_file = "%s-%s%s" % (base, region, ext)
    chrom, start, end = region.split("-")
    start = int(start) - 1
    end = int(end)
    if not os.path.exists(out_file):
        index = pyfaidx.Fasta(in_file)
        with open(out_file, "w") as out_handle:
            for recid in index.keys():
                out_handle.write(repr(index[recid][start:end]) + "\n")
    return out_file

def _subset_ref_file(ref_file, out_dir, region):
    chrom, start, end = region.split("-")
    base, ext = os.path.splitext(os.path.basename(ref_file))
    out_file = os.path.join(out_dir, "%s-%s%s" % (base, region, ext))
    region = "%s:%s-%s" % (chrom, int(start), int(end))
    cmd = "samtools faidx {ref_file} {region} | sed 's/{region}/{chrom}/' > {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def _index_ref_file(ref_file):
    out_file = "%s.fai" % ref_file
    if not os.path.exists(out_file):
        cmd = ["samtools", "faidx", ref_file]
        subprocess.check_call(cmd)
    bwa_index = "%s.bwt" % ref_file
    if not os.path.exists(bwa_index):
        cmd = ["bwa", "index", ref_file]
        subprocess.check_call(cmd)
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
