#!/usr/bin/env python
"""Prepare input configuration files for alignment and QC of reads.
"""
import csv
import glob
import os
import re
import subprocess
import sys

def main(template_file, *input_dirs):
    name_re = re.compile("-6\d\d\d-")
    fnames = []
    for input_dir in input_dirs:
        fnames.extend([f for f in glob.glob(os.path.join(input_dir, "*.fastq.gz"))
                       if not os.path.basename(f).startswith("Undetermined")])
    out_file = os.path.join(os.path.dirname(template_file), "li_hiv_call3.csv")

    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["samplename", "description"])
        for fname in fnames:
            descr = name_re.split(os.path.basename(fname))[0]
            writer.writerow([os.path.basename(fname), descr])
    cmd = ["bcbio_nextgen.py", "-w", "template", template_file, out_file] + fnames
    subprocess.check_call(cmd)

if __name__ == "__main__":
    main(*sys.argv[1:])
