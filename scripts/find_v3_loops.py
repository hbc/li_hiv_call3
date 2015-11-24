#!/usr/bin/env python
"""Find V3 loops from a group of targeted sequencing samples.
"""
import collections
import os
import subprocess
import sys

from Bio import SeqIO

def main(ref_v3_file, *fastq_files):
    v3_seqs = _read_v3_file(ref_v3_file)
    for base, fq1, fq2 in _pair_fastq_files(fastq_files):
        merge_fq = _merge_fastq(fq1, fq2, base)
        counts = _count_v3_seqs(merge_fq, v3_seqs)
        print base, counts

def _count_v3_seqs(merge_fq, v3_seqs):
    counts = collections.defaultdict(int)
    for rec in SeqIO.parse(merge_fq, "fastq"):
        for v3_name, v3_seq in v3_seqs.items():
            if str(rec.seq).find(v3_seq[:50]) >= 0 or str(rec.seq).find(v3_seq[-50:]) >= 0:
                counts[v3_name] += 1
    return dict(counts)

def _pair_fastq_files(in_files):
    out = collections.defaultdict(list)
    for f in in_files:
        base = os.path.basename(f).split("_")[0]
        out[base].append(f)
    for base, fnames in out.items():
        if len(fnames) == 2:
            fnames.sort()
            yield base, fnames[0], fnames[1]

def _merge_fastq(fq1, fq2, base):
    """Extract paired end reads and merge using pear.
    """
    out_dir = "merged"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_base = os.path.join(out_dir, "%s-merge" % (base))
    out_file = out_base + ".assembled.fastq"
    if not os.path.exists(out_file):
        cmd = "pear -q 3 -v 40 -t 75 -n 90 -u 1 -f {fq1} -r {fq2} -o {out_base}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def _read_v3_file(in_file):
    out = {}
    for rec in SeqIO.parse(in_file, "fasta"):
        out[rec.id] = str(rec.seq)
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])