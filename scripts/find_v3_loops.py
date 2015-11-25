#!/usr/bin/env python
"""Find V3 loops from a group of targeted sequencing samples.

Jon and Zixin;
Thanks for the discussion and sorry for not being more clear in the process
used. You're right on in that these are exact matches
"""
import collections
import contextlib
import os
import subprocess
import sys

import pysam

from Bio import pairwise2, SeqIO

def main(ref_v3_file, *fastq_files):
    v3_seqs = _read_v3_file(ref_v3_file)
    for base, fq1, fq2 in _pair_fastq_files(fastq_files):
        merge_fq = _merge_fastq(fq1, fq2, base)
        #counts = _count_v3_seqs(merge_fq, v3_seqs)
        counts = _align_to_v3_seqs(merge_fq, ref_v3_file)
        print base, counts

def _align_to_v3_seqs(merge_fq, ref_file):
    out_file = "%s.bam" % os.path.splitext(merge_fq)[0]
    counts = collections.defaultdict(int)
    if not os.path.exists(out_file):
        cmd = "bwa index {ref_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
        cmd = "bwa mem {ref_file} {merge_fq} | samtools sort -T {out_file}-tmp -O bam -o {out_file} /dev/stdin"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    with contextlib.closing(pysam.Samfile(out_file)) as work_bam:
        for read in work_bam:
            if read.reference_id >= 0:
                counts[read.reference_name] += 1
    return dict(counts)

def _count_v3_seqs(merge_fq, v3_seqs):
    counts = collections.defaultdict(int)
    for rec in SeqIO.parse(merge_fq, "fastq"):
        scores = []
        for v3_name, v3_seq in v3_seqs.items():
            if str(rec.seq).find(v3_seq[:50]) >= 0 or str(rec.seq).find(v3_seq[-50:]) >= 0:
                counts[v3_name] += 1
            #score = _align_two(str(rec.seq), v3_seq)
            #scores.append((score, v3_name))
        #scores.sort(reverse=True)
        #best_score, best_name = scores[0]
        #if best_score > 50:
        #    counts[best_name] += 1
    return dict(counts)

def _align_two(s1, s2):
    alignments = pairwise2.align.localms(s1, s2, 2, -1, -1, -0.5, one_alignment_only=True)
    matches = 0
    if len(alignments) > 0:
        for b1, b2 in zip(alignments[0][0], alignments[0][1]):
            if b1 != "-" and b1 == b2:
                matches += 1
    return matches

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