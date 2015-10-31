#!/usr/bin/env python
"""Convert TherapyEdge outputs to VCF to do comparisons against the truth set.
"""
import os
import sys

import pandas

def main(in_file):
    df = pandas.read_csv(in_file, sep=";")
    df = df[df["Passed Expert System"] == "Yes"]
    gdfs = {}
    for region, gdf in df.groupby("Protein"):
        gdfs[region] = gdf

    out_file = "%s.vcf" % os.path.splitext(in_file)[0]
    with open(out_file, "w") as out_handle:
        out_handle.write("##fileformat=VCFv4.1\n")
        out_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcall\n")
        print gdfs.keys()
        for region, start in [("PROT", 2252), ("RT", 2549), ("INT", 4229), ("GP120", 6224)]:
            print gdfs[region].head()
            for pos, posdf in gdfs[region].groupby("Position"):
                alts = []
                freqs = []
                for i, row in posdf.iterrows():
                    ref, alt = row["Mutation"].split(str(row["Position"]))
                    alts.append(alt)
                    freqs.append(str(row["Prevalence"]))
                line = ["HIVHXB2CG", start + pos, ".", ref, ",".join(alts), ".", ".",
                        "AF=%s;FREQ=%s" % (max([float(f) for f in freqs]) / 100.0, ",".join(freqs)),
                        "GT", "0/1"]
                out_handle.write("\t".join(str(x) for x in line) + "\n")

if __name__ == "__main__":
    main(*sys.argv[1:])