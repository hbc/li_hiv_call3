#!/usr/bin/env python
"""Plot distribution of allele frequencies.

Usage:
  plot_af_distribution.py region.bed sample1.vcf sample2.vcf ...
"""
import os
import sys

import vcf

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(region_file, *vcf_files):
    regions = read_region_file(region_file)
    sample_afs = []
    for vcf_file in vcf_files:
        sample_afs.extend(parse_vcf(vcf_file, regions))
    df = pd.DataFrame(sample_afs)
    df = df[(df.af > 1.0) & (df.af < 25.0)]
    sns.set(style="dark")
    g = sns.FacetGrid(df, row="region", col="sample", sharey=False, margin_titles=True)
    g.map(plt.hist, "af", bins=50)
    #g.map(sns.distplot, "af")
    g.set(xlim=(0.0, 25.0))
    g.fig.set_size_inches(8, 11)
    g.fig.savefig("seracare_feb2016_freqs.pdf")

def parse_vcf(in_file, regions):
    out = []
    sample = os.path.basename(in_file).split("-")[0]
    with open(in_file) as in_handle:
        reader = vcf.Reader(in_handle)
        for rec in reader:
            if len(rec.FILTER) == 0:
                region = None
                for chrom, start, end, name in regions:
                    if chrom == rec.CHROM and rec.POS >= start and rec.POS <= end:
                        region = name
                        break
                out.append({"sample": sample, "region": region, "af": rec.INFO["AF"] * 100.0})
    return out

def read_region_file(in_file):
    regions = []
    with open(in_file) as in_handle:
        for line in in_handle:
            chrom, start, end, name = line.strip().split("\t")
            regions.append((chrom, int(start), int(end), name))
    return regions

if __name__ == "__main__":
    main(*sys.argv[1:])
