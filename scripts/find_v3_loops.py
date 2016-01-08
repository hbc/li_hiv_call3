#!/usr/bin/env python
"""Find V3 loops from a group of targeted sequencing samples.
"""
import collections
import contextlib
import itertools
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
        #print "counting", base, counts
        counts = _align_to_v3_seqs(merge_fq, ref_v3_file)
        print "aligning", base, counts
        _extract_seqs_from_align(merge_fq, base, "%s-translate.fa" % os.path.splitext(ref_v3_file)[0])

def _extract_seqs_from_align(merge_fq, base, ref_v3_p_file):
    bam_file = "%s.bam" % os.path.splitext(merge_fq)[0]
    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    with contextlib.closing(pysam.Samfile(bam_file)) as work_bam:
        for read in work_bam:
            if read.reference_id >= 0:
                if len(read.query_alignment_sequence) > 70:
                    counts[read.reference_name][read.query_alignment_sequence] += 1
    match_dir = "matches"
    if not os.path.exists(match_dir):
        os.makedirs(match_dir)
    for match_name in counts.keys():
        total = 0
        out_file = os.path.join(match_dir, "%s_match_%s.fa" % (base, match_name))
        with open(out_file, "w") as out_handle:
            for i, (count, seq) in enumerate(sorted([[c, seq] for seq, c in counts[match_name].items()], reverse=True)):
                total += count
                out_handle.write(">%s-count=%s\n%s\n" % (i, count, seq))
        print base, match_name, total, len(counts[match_name])
        _translate_seq(out_file, base, match_name, ref_v3_p_file)

def _translate_seq(in_file, base, match_name, ref_v3_p_file):
    v3_prots = _read_v3_file(ref_v3_p_file)
    prot_file = "%s-orfs.fa" % os.path.splitext(in_file)[0]
    out_file = "%s-protein.fa" % os.path.splitext(in_file)[0]
    cmd = "getorf -minsize 20 -noreverse -outseq {prot_file} -sequence {in_file}"
    print cmd.format(**locals())
    subprocess.check_call(cmd.format(**locals()), shell=True)
    with open(out_file, "w") as out_handle:
        SeqIO.write(_find_best_aa(SeqIO.parse(prot_file, "fasta"), v3_prots[match_name]), out_handle, "fasta")

def _find_best_aa(recs, v3_seq):
    for cur_name, rec_group in itertools.groupby(recs, lambda x: x.id.split("_")[0]):
        max_score = 0
        best = None
        for rec in rec_group:
            if _align_two(str(rec.seq), v3_seq) > max_score:
                best = rec
        best.id = cur_name
        best.description = ""
        yield best

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
            if str(rec.seq).find(v3_seq[20:45]) >= 0 and str(rec.seq).find(v3_seq[-45:-20]) >= 0:
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