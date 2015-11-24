#!/usr/bin/env python
"""Collapse reads on V3 loop and organize for calling.
"""
import os
import subprocess
import sys

import yaml

def main(config_file, bam_file):
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    #v3_region = _read_regions_file(config["regions"])
    v3_region = "HIVHXB2CG:7110-7217"
    [fq1, fq2], [full_fq1, full_fq2] = _select_fastq_in_region(bam_file, v3_region)
    merged_fq = _merge_fastq(fq1, fq2, config)
    print calculate_kmers(merged_fq)
    print align_to_control(full_fq1, full_fq2, config)

def align_to_control(fq1, fq2, config):
    """Align fastqs to control file, counting hits per input.
    """
    ref_file = config["controls"]["V3-Loop"]
    if not os.path.exists(ref_file + ".bwt"):
        cmd = "bwa index {ref_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    out_file = "%s-v3_counts.txt" % os.path.splitext(fq1)[0]
    if not os.path.exists(out_file):
        cmd = ("bwa mem {ref_file} {fq1} {fq2} | samtools view -F 4 - | cut -f 3 | sort | uniq -c > {out_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def calculate_kmers(merged_fq):
    """Use jellyfish to count kmers in the input file.
    """
    jf_file = "%s.jf" % os.path.splitext(merged_fq)[0]
    if not os.path.exists(jf_file):
        cmd = "jellyfish count -s 100M -m 105 -o {jf_file} {merged_fq}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    kmer_file = "%s.kmer" % os.path.splitext(merged_fq)[0]
    if not os.path.exists(kmer_file):
        cmd = "jellyfish dump -c -t -L 100 {jf_file} | sort -k 2nr > {kmer_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return kmer_file

def _read_regions_file(in_file):
    with open(in_file) as in_handle:
        for line in in_handle:
            contig, start, end, name = line.strip().split()
            if name.lower().startswith("v3"):
                return "%s:%s-%s" % (contig, start, end)

def _merge_fastq(fq1, fq2, config):
    """Extract paired end reads and merge using pear.
    """
    base_out = os.path.splitext(fq1)[0].replace("-1", "-pear")
    merged_out = "%s.assembled.fastq" % base_out
    if not os.path.exists(merged_out):
        cmd = ["pear", "-q", str(config["params"]["pear"]["quality_thresh"]),
               "-v", str(config["params"]["pear"]["min_overlap"]),
               "-t", str(config["params"]["pear"]["min_trim_length"]),
               "-u", str(config["params"]["pear"]["max_uncalled_base"]),
               "-f", fq1, "-r", fq2, "-o", base_out]
        subprocess.check_call(cmd)
    return merged_out

def _select_fastq_in_region(bam_file, region):
    """Extract paired end reads where either read overlaps the region of interest.
    """
    work_dir = os.path.join(os.getcwd(), "v3_fastq")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    full_fq1, full_fq2 = _select_full_fastqs(bam_file, work_dir)
    base, ext = os.path.splitext(os.path.basename(bam_file))
    base = os.path.join(work_dir, base)

    region_str = region.replace(":", "_").replace("-", "_")
    fq1 = "%s-%s-1.fastq" % (base, region_str)
    fq2 = "%s-%s-2.fastq" % (base, region_str)
    name_file = "%s-%s-names.txt" % (base, region_str)
    keep_file = "%s-keep%s" % os.path.splitext(name_file)
    if not os.path.exists(keep_file):
        cmd = ("samtools view {bam_file} {region} | cut -f 1 | sort | uniq > {name_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
        with open(name_file) as in_handle:
            with open(keep_file, "w") as out_handle:
                for line in in_handle:
                    out_handle.write("%s/1\n" % line.strip())
                    out_handle.write("%s/2\n" % line.strip())
    for orig_file, out_file in [(full_fq1, fq1), (full_fq2, fq2)]:
        if not os.path.exists(out_file):
            cmd = "seqtk subseq {orig_file} {keep_file} > {out_file}"
            subprocess.check_call(cmd.format(**locals()), shell=True)
    return [fq1, fq2], [full_fq1, full_fq2]

def _select_full_fastqs(bam_file, work_dir):
    base, ext = os.path.splitext(os.path.basename(bam_file))
    base = os.path.join(work_dir, base)
    full_fq1 = "%s-1.fastq" % base
    full_fq2 = "%s-2.fastq" % base
    if not os.path.exists(full_fq1):
        cmd = ("bamtofastq F={full_fq1} F2={full_fq2} S=/dev/null O=/dev/null O2=/dev/null collate=1 "
               "filename={bam_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return full_fq1, full_fq2

if __name__ == "__main__":
    main(*sys.argv[1:])