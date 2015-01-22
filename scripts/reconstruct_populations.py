#!/usr/bin/env python
"""Reconstruct HIV populations from aligned BAM files in regions of interest.

Requires:
  bedtools
  ViQuaS

Usage:
    reconstruct_populations.py <config_file> [<BAM files>]
"""
import os
import sys
import subprocess

import yaml

def main(config_file, *bam_files):
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    regions = _subset_region_file(config["regions"])
    for bam_file in bam_files:
        for region, region_file in regions:
            subbam_file = _select_regions(bam_file, region_file, region)
            recon_dir = _run_viquas(subbam_file, config["ref_file"])
            print recon_dir

def _run_viquas(bam_file, ref_file):
    out_dir = "%s-viquas" % os.path.splitext(bam_file)[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    orig_dir = os.getcwd()
    os.chdir(out_dir)

    cmd = "Rscript ViQuaS.R {ref_file} {bam_file} 0 1 1 {size}"
    os.chdir(orig_dir)

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
