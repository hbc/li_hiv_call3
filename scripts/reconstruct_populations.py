#!/usr/bin/env python
"""Reconstruct HIV populations from aligned BAM files in regions of interest.

Requires:
  bedtools
  samtools
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

import yaml
import joblib

def main(cores, config_file, *bam_files):
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    regions = _subset_region_file(config["regions"])
    to_run = []
    for bam_file in bam_files:
        for region, region_file in regions:
            subbam_file = _select_regions(bam_file, region_file, region)
            to_run.append((subbam_file, config, region))
    for recon_file in joblib.Parallel(int(cores))(joblib.delayed(_run_viquas)(f, c, r) for f, c, r in to_run):
        print recon_file

def _run_viquas(bam_file, config, region):
    out_dir = "%s-viquas" % os.path.splitext(bam_file)[0]
    out_file = os.path.join(out_dir, "ViQuaS-Spectrum.fa")
    bam_file = os.path.join(os.getcwd(), bam_file)
    ref_file = os.path.join(os.getcwd(), config["ref_file"])
    viquas_dir = os.path.join(os.getcwd(), config["viquas_dir"])
    # parameter settings for Illumina data from ViQuaS paper
    o = 5
    r = 0.7
    chrom, start, end = region.split("-")
    size = int(end) - int(start)
    if not os.path.exists(out_file):
        with _chdir(out_dir):
            with _copy_viquas(viquas_dir):
                cmd = "Rscript ViQuaS.R {ref_file} {bam_file} {o} '' 1 {size}"
                subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

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

if __name__ == "__main__":
    main(*sys.argv[1:])
