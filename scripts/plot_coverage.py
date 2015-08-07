#!/usr/bin/env python
"""Plot coverage before and after normalization for control samples.
"""
import os
import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(stats_file):
    out_file = "%s.pdf" % os.path.splitext(stats_file)[0]
    data = {"bam": [], "start": [], "coverage": []}
    with open(stats_file) as in_handle:
        cur_bam = None
        for line in in_handle:
            if line.find(".bam") >= 0:
                cur_bam = line.strip().replace("-ready", "").replace(".bam", "")
            elif not line.startswith("#") and cur_bam:
                _, start, _, count, coverage, sample = line.strip().split()[:6]
                data["bam"].append(cur_bam)
                data["start"].append(int(start))
                data["coverage"].append(float(coverage))
    df = pd.DataFrame(data)
    g = sns.FacetGrid(df, col="bam", col_wrap=2)
    g.map(plt.plot, "start", "coverage")
    g.fig.set_size_inches(10, 8)
    g.fig.savefig(out_file)

if __name__ == "__main__":
    main(*sys.argv[1:])
