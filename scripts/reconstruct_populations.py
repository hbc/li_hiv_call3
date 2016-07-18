#!/usr/bin/env python
"""Reconstruct HIV populations from aligned BAM files in regions of interest.

Requires:
  bedtools
  samtools
  bfc (https://github.com/lh3/bfc)
  pear (http://www.exelixis-lab.org/web/software/pear)
  bwa
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
from __future__ import print_function
import collections
import contextlib
import itertools
import os
import pprint
import re
import shutil
import sys
import subprocess

from Bio import pairwise2, SeqIO
import joblib
import numpy
import pyfaidx
import pysam
import yaml

def main(cores, config_file, *bam_files):
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    regions = _subset_region_file(config["regions"])
    to_prep = []
    to_call = []
    if config["params"]["caller"] == "lofreq":
        if config.get("control_vcf"):
            control_bam = bam_files[0]
            for bam_file in bam_files[1:]:
                to_call.append((bam_file, control_bam, config["ref_file"], config["regions"], config))
        else:
            for bam_file in bam_files:
                to_call.append((bam_file, None, config["ref_file"], config["regions"], config))
    else:
        for bam_file in bam_files:
            to_call.append((bam_file, config["ref_file"], config["regions"], config))
            for region, name, region_file in regions:
                to_prep.append((bam_file, region, name, region_file, config))

    if config["params"]["caller"] == "freebayes":
        called = joblib.Parallel(int(cores))(joblib.delayed(_run_freebayes)(*args) for args in to_call)
    elif config["params"]["caller"] == "lofreq":
        called = joblib.Parallel(int(cores))(joblib.delayed(_run_lofreq)(*args) for args in to_call)
    elif config["params"]["caller"] == "vardict":
        called = joblib.Parallel(int(cores))(joblib.delayed(_run_vardict)(*args) for args in to_call)
    elif config["params"]["caller"] == "viquas":
        to_run = joblib.Parallel(int(cores))(joblib.delayed(_run_prep)(*args) for args in to_prep)
        def by_size((f, c, r, n)):
            return os.path.getsize(f)
        to_run.sort(key=by_size)
        for recon_file in joblib.Parallel(int(cores))(joblib.delayed(_run_viquas)(f, c, r, n) for f, c, r, n in to_run):
            print(recon_file)
        called = None
    if called:
        for vcf_file in called:
            if "Control" in vcf_file:
                _summarize_calls(vcf_file, config["ref_file"], config["control_vcf"],
                                 config["val_regions"], config["params"]["validation"])

def _run_prep(bam_file, region, name, region_file, config):
    prep_fq1, prep_fq2 = _select_fastq_in_region(bam_file, region, region_file)
    merged_bam_file = _run_cleaning(prep_fq1, prep_fq2, region_file, config)
    _plot_coverage(merged_bam_file)
    return (merged_bam_file, config, region, name)

def _run_cleaning(prep_fq1, prep_fq2, region_file, config):
    fai_file = _index_ref_file(config["ref_file"])
    if "bfc" in config["params"]["prep"]:
        corr_fq1 = _correct_fastq(prep_fq1, fai_file, config)
        corr_fq2 = _correct_fastq(prep_fq2, fai_file, config)
        prep_fq1, prep_fq2 = _ensure_pairs(corr_fq1, corr_fq2)
    if "pear" in config["params"]["prep"]:
        # _check_pair_overlap(prep_fq1, prep_fq2, config["ref_file"], config)
        prep_fq1 = _merge_fastq(prep_fq1, prep_fq2, config)
        prep_fq2 = None
    return _realign_merged(prep_fq1, prep_fq2, region_file, config["ref_file"])

def _run_lofreq(bam_file, control_bam, ref_file, region_file, config):
    """Perform pooled variant calling with lofreq.
    """
    params = config["params"]["lofreq"]
    work_dir = "lofreq"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if params["prep"].startswith("gatk"):
        if params["prep"] == "gatk_viterbi":
            realign_bam_file = os.path.join(work_dir, "%s-realign%s" % os.path.splitext(os.path.basename(bam_file)))
            cmd = ("lofreq viterbi -f {ref_file} {bam_file} | "
                   "samtools sort -O bam -T {realign_bam_file}-tmp -o {realign_bam_file}")
            if not os.path.exists(realign_bam_file):
                subprocess.check_call(cmd.format(**locals()), shell=True)
            bai_file = realign_bam_file + ".bai"
            if not os.path.exists(bai_file):
                cmd = "samtools index {realign_bam_file}"
                subprocess.check_call(cmd.format(**locals()), shell=True)
        else:
            realign_bam_file = bam_file
        recal_file = _get_bqsr_recal_file_control(control_bam, work_dir, ref_file, region_file,
                                                  config["control_vcf"])
        prep_bam_file = os.path.join(work_dir, "%s-prep%s" % os.path.splitext(os.path.basename(realign_bam_file)))
        cmd = ("java -jar GenomeAnalysisTK.jar -T PrintReads -R {ref_file} "
               "-I {realign_bam_file} -L {region_file} -BQSR {recal_file} -o {prep_bam_file}")
        if not os.path.exists(prep_bam_file):
            subprocess.check_call(cmd.format(**locals()), shell=True)
        extra_args = "--call-indels"
    elif params["prep"] == "lofreq":
        prep_bam_file = os.path.join(work_dir, "%s-prep%s" % os.path.splitext(os.path.basename(bam_file)))
        cmd = ("lofreq viterbi -f {ref_file} {bam_file} | "
               "lofreq indelqual --dindel -f {ref_file} - | "
               "samtools sort -O bam -T {prep_bam_file}-tmp -o {prep_bam_file}")
        if not os.path.exists(prep_bam_file):
            subprocess.check_call(cmd.format(**locals()), shell=True)
        extra_args = "--call-indels"
    else:
        prep_bam_file = bam_file
        extra_args = ""
    bai_file = prep_bam_file + ".bai"
    if not os.path.exists(bai_file):
        cmd = "samtools index {prep_bam_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    out_vcf = os.path.join(work_dir, "%s.vcf" % os.path.splitext(os.path.basename(prep_bam_file))[0])
    if not os.path.exists(out_vcf) or os.path.getsize(out_vcf) == 0:
        sb_thresh = params["sb_thresh"]
        # cmd = ("lofreq call-parallel --pp-threads 4 {prep_bam_file} -f {ref_file} -l {region_file} -N {extra_args} |"
        #        "lofreq filter -B {sb_thresh} --sb-incl-indels | vt normalize -r {ref_file} - > {out_vcf}")
        downsample_n = config["params"]["lofreq"].get("downsample")
        downsample = "-d %s" % (downsample_n) if downsample_n else ""
        cmd = ("lofreq call-parallel --pp-threads 8 {prep_bam_file} {downsample} -f {ref_file} -l {region_file} {extra_args} "
               "--no-default-filter > {out_vcf}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return _filter_lofreq(out_vcf)

def _filter_lofreq(out_vcf):
    filter_vcf = "%s-filter%s" % os.path.splitext(out_vcf)
    if not os.path.exists(filter_vcf) or os.path.getsize(filter_vcf) == 0:
        cmd = ("bcftools filter -o {filter_vcf} -s LowQual "
               "-e 'QUAL < 200 & AF < 0.01' "
               "{out_vcf}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return filter_vcf

def _get_bqsr_recal_file_control(bam_file, work_dir, ref_file, region_file, known_vcf):
    """Retrieve a BQSR recalibration file using control BAM file with a truth VCF.

    Trains recalibration against the actual low frequency variants in the input.
    """
    recal_file = os.path.join(work_dir, "%s-prep.grp" % os.path.splitext(os.path.basename(bam_file))[0])
    cmd = ("java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R {ref_file} "
            "-I {bam_file} -L {region_file} -knownSites {known_vcf} -o {recal_file}")
    if not os.path.exists(recal_file):
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return recal_file

def _get_bqsr_recal_file_naive(bam_file, work_dir, ref_file, region_file):
    """Naive approach to get a recalibration file using samtools mpileup calls as truth set.
    """
    known_vcf = os.path.join(work_dir, "%s-prep-known.vcf.gz" % os.path.splitext(os.path.basename(bam_file))[0])
    cmd = ("samtools mpileup -uf {ref_file} {bam_file} | bcftools call -mv -Oz > {known_vcf}")
    if not os.path.exists(known_vcf):
        subprocess.check_call(cmd.format(**locals()), shell=True)
    if not os.path.exists(known_vcf + ".tbi"):
        cmd = "tabix -p vcf -f {known_vcf}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return _get_bqsr_recal_file_control(bam_file, work_dir, ref_file, region_file, known_vcf)

def _run_freebayes(bam_file, ref_file, region_file, config):
    """Merge and call haplotype variants with FreeBayes.
    """
    work_dir = "freebayes"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if config["params"]["freebayes"]["clean_bam"]:
        full_fq1, full_fq2 = _select_full_fastqs(bam_file, work_dir)
        full_bam = _run_cleaning(full_fq1, full_fq2, None, config)
    else:
        full_bam = bam_file
    vcf_file = _freebayes_call(full_bam, ref_file, region_file, config)
    return vcf_file

def _run_viquas(bam_file, config, region, region_name):
    out_dir = os.path.join(os.getcwd(), "%s-viquas" % os.path.splitext(bam_file)[0])
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file = os.path.join(out_dir, "ViQuaS-Spectrum.fa")
    bam_file = os.path.join(os.getcwd(), bam_file)
    viquas_dir = os.path.join(os.getcwd(), config["viquas_dir"])
    ref_file = os.path.join(os.getcwd(), config["ref_file"])
    o = config["params"]["viquas"]["o"]
    r = config["params"]["viquas"]["r"]
    chrom, start, end = region.split("-")
    size = int(end) - int(start)
    if not os.path.exists(out_file):
        with _chdir(out_dir):
            with _copy_viquas(viquas_dir):
                cmd = "Rscript ViQuaS.R {ref_file} {bam_file} {o} {r} 0 {size}"
                subprocess.check_call(cmd.format(**locals()), shell=True)
    called_file = _subset_viquas_file(out_file, region, config["ref_file"])
    if _is_control(called_file):
        _evaluate_control(called_file, config["controls"][re.split("_\d+", region_name)[0]])
    return called_file

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

def _fix_r_parsing(line):
    """Fix parsing of the 'r' parameter to ensure it handles floats.
    """
    if line.startswith("r = as."):
        return "r = as.numeric(ip[4])\n"
    else:
        return line

@contextlib.contextmanager
def _copy_viquas(viquas_dir):
    """ViQuaS requires running within the install directory.
    This copies it into the current directory so multiple instances can be
    run simultaneously.
    """
    patch_fns = {"ViQuaS.R": _fix_r_parsing}
    files = ["ViQuaS.R"]
    dirs = ["tools", "viquas_files"]
    for fname in files:
        orig_file = os.path.join(viquas_dir, fname)
        out_file = os.path.join(os.getcwd(), fname)
        if not os.path.exists(out_file):
            with open(out_file, "w") as out_handle:
                with open(orig_file) as in_handle:
                    for line in in_handle:
                        if fname in patch_fns:
                            line = _fix_r_parsing(line)
                        out_handle.write(line)
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

def _correct_fastq(in_file, fai_file, config):
    """Error correct input fastq file with bfc.
    """
    gsize = 0
    with open(fai_file) as in_handle:
        for line in in_handle:
            gsize += int(line.split()[1])
    out_file = "%s-corrected%s" % os.path.splitext(in_file)
    min_kmer_cov = config["params"]["bfc"]["min_kmer_cov"]
    if not os.path.exists(out_file):
        cmd = "bfc -s {gsize} -E -1 -c {min_kmer_cov} {in_file} > {out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def _ensure_pairs(f1, f2):
    """Ensure only pairs are retained with both ends.
    """
    out_files = []
    for filter_f, cmp_f, new_ext, orig_ext in [(f1, f2, "/1", "/2"), (f2, f1, "/2", "/1")]:
        out_file = "%s-wboth%s" % os.path.splitext(filter_f)
        if not os.path.exists(out_file):
            keep_file = "%s-names.txt" % os.path.splitext(out_file)[0]
            with open(cmp_f) as in_handle:
                with open(keep_file, "w") as out_handle:
                    for line in in_handle:
                        if line.startswith("@"):
                            out_handle.write(line[1:].replace(orig_ext, new_ext))
            cmd = "seqtk subseq {filter_f} {keep_file} > {out_file}"
            subprocess.check_call(cmd.format(**locals()), shell=True)
        out_files.append(out_file)
    return out_files

# -- Validation

def _summarize_calls(vcf_file, ref_file, control_file, region_file, params):
    """Summarize calls against standard validation VCF, not requiring phasing.
    """
    want_pos = _read_regions(region_file, params["end_buffer"])
    control_calls = _read_control(control_file, want_pos)

    tps = []
    tp_freqs = collections.defaultdict(int)
    fps = []
    fps_lowfreq = []
    abovefreq = {"TP": 0, "FP": 0, "FN": 0}
    sb_tps = 0
    sb_fps = 0
    with pysam.VariantFile(vcf_file) as in_bcf:
        for rec in in_bcf:
            if rec.pos in want_pos[rec.chrom] and "PASS" in rec.filter:
                for alt in rec.alts:
                    key = (rec.chrom, rec.pos, rec.ref, alt)
                    freq = rec.info["AF"] * 100.0
                    if key in control_calls:
                        tp_freqs[control_calls[key]] += 1
                        tps.append(control_calls[key])
                        freq = control_calls[key]
                        del control_calls[key]
                        if rec.info.get("SB", 0) > params["sb_thresh"]:
                            print(dict(rec.info))
                            sb_tps += 1
                        if freq >= params["fp_thresh"]:
                            abovefreq["TP"] += 1
                    elif rec.info.get("SB", 0) > params["sb_thresh"]:
                        sb_fps += 1
                    elif float(rec.info["AF"]) >= params["fp_thresh"] / 100.0:
                        fps.append((key, dict(rec.info)))
                        abovefreq["FP"] += 1
                    else:
                        fps_lowfreq.append((key, dict(rec.info)))
    for freq in control_calls.itervalues():
        if freq >= params["fp_thresh"]:
            abovefreq["FN"] += 1
    print("Strand bias filter", sb_tps, sb_fps)
    print()
    print("TP", len(tps), "\nFP", len(fps) + len(fps_lowfreq), "\nFN", len(control_calls),
          "\nFP lowfreq (< %.2f)" % params["fp_thresh"], len(fps_lowfreq))
    print()
    for key, count in abovefreq.items():
        print("Above frequency (%.2f) %s: %s" % (params["fp_thresh"], key, count))
    print()
    fn_freqs = collections.defaultdict(int)
    print("High frequency FNs (> %.2f)" % params["fn_highfreq_thresh"])
    for key, freq in sorted(control_calls.items(), key=lambda x: x[0][1]):
        fn_freqs[freq] += 1
        if freq > params["fn_highfreq_thresh"]:
            print(key, freq)
    print("TP/FN frequencies")
    freqs = list(sorted(set(tp_freqs.keys() + fn_freqs.keys()), reverse=True))
    fp_freqs = collections.defaultdict(int)
    for fp in sorted(fps_lowfreq + fps):
        cur_freq = float(fp[1]["AF"]) * 100.0
        freq = min(freqs, key=lambda x: abs(x - cur_freq))
        fp_freqs[freq] += 1
    print("| freq | TP | FN | FP |")
    for freq in freqs:
        print("| %.1f | %s | %s | %s |" % (freq, tp_freqs[freq], fn_freqs[freq], fp_freqs[freq]))
    print("All false positives and negatives")
    wrong = []
    for fp in sorted(fps_lowfreq + fps):
        wrong.append((fp[0], "%.1f" % (float(fp[1]["AF"]) * 100.0), "FP"))
    for fn in control_calls.items():
        wrong.append((fn[0], "%.1f" % fn[1], "FN"))
    for x in sorted(wrong):
        if x[0][1] in want_pos[x[0][0]]:
            print(x)

def _read_regions(region_file, end_buffer):
    """Read supported positions by contig, providing an end buffer around PCR amplification regions.
    """
    regions = collections.defaultdict(set)
    with open(region_file) as in_handle:
        for line in in_handle:
            chrom, start, end = line.strip().split()[:3]
            regions[chrom] |= set(range(int(start) + end_buffer, int(end) - end_buffer))
    return regions

def _read_control(control_file, want_pos):
    """Read control VCF into a dictionary to use for lookups.
    """
    out = {}
    with pysam.VariantFile(control_file) as in_vcf:
        for rec in in_vcf:
            if rec.pos in want_pos[rec.chrom]:
                freqs = rec.info["FREQ"]
                if not isinstance(freqs, (list, tuple)):
                    freqs = [freqs]
                assert len(freqs) == len(rec.alts), (rec.chrom, rec.pos, freqs, rec.alts)
                for alt, freq, in zip(rec.alts, freqs):
                    key = (rec.chrom, rec.pos, rec.ref, alt)
                    out[key] = freq
    return dict(out)

def _summarize_calls_haplotype(vcf_file, ref_file, region_file, config):
    """Summarize calls from a VCF file compared against known controls, requiring long haplotypes.
    """
    sizes = []
    cmanage = ControlManager(region_file, ref_file, config)
    with pysam.VariantFile(vcf_file) as in_vcf:
        for rec in in_vcf:
            counts = [int(x) for x in
                      (rec.info["AO"] if isinstance(rec.info["AO"], (list, tuple)) else [rec.info["AO"]])]
            freqs = [float(x) / float(rec.info["DP"]) * 100.0 for x in counts]
            allele_freqs = [(a, f, c) for f, a, c in sorted(zip(freqs, rec.alts, counts), reverse=True)]
            if len(allele_freqs[0][0]) > 10:
                sizes.append(len(allele_freqs[0][0]))
            cmanage = _check_haplotype_matches(cmanage, rec.chrom, rec.start, rec.ref, allele_freqs,
                                               ref_file, region_file, config)
    print("Wrong frequency")
    pprint.pprint(dict(cmanage.evals["wrongfreq"]))
    print("True positives")
    _print_eval_summary(dict(cmanage.evals["tps"]), 50.0)
    pprint.pprint(dict(cmanage.evals["tps"]))
    print("False positives")
    _print_eval_summary(dict(cmanage.evals["fps"]), 0.4)
    pprint.pprint(dict(cmanage.evals["fps"]))
    print("False negatives")
    for control in cmanage.get_all_controls():
        print(control.name, control.get_no_matches())
    if len(sizes) > 0:
        print("Sizes")
        print(numpy.median(sizes), max(sizes), min(sizes))

def _print_eval_summary(vals, thresh):
    out = []
    above = 0
    for freq, count in sorted(vals.items(), key=lambda x: float(x[0])):
        if float(freq) > thresh:
            above += 1
        else:
            out.append("%s: %s" % (freq, count))
    out.append("above %s: %s" % (thresh, above))
    print("\n".join(out) + "\n")

def _get_control_file(in_chrom, in_pos, region_file, config):
    with open(region_file) as in_handle:
        for line in in_handle:
            chrom, start, end, name = line.rstrip().split("\t")
            if in_chrom == chrom and in_pos >= int(start) and in_pos < int(end):
                return name, (chrom, int(start), int(end)), config["controls"][name]
    raise ValueError("Did not find control for %s %s" % (in_chrom, in_pos))

class ControlManager:
    def __init__(self, region_file, ref_file, config):
        self._region_file = region_file
        self._ref_file = ref_file
        self._config = config
        self._controls = {}
        self.evals = {"tps": collections.defaultdict(int),
                      "fps": collections.defaultdict(int),
                      "wrongfreq": collections.defaultdict(int)}

    def get_all_controls(self):
        out = []
        for control_file, controls in self._controls.items():
            out.extend(controls)
        return out

    def get_controls(self, chrom, start):
        control_name, control_coords, control_file = _get_control_file(chrom, start, self._region_file,
                                                                       self._config)
        if control_file in self._controls:
            return self._controls[control_file]
        else:
            controls = _prep_control_alignments(control_file, control_coords, self._ref_file)
            self._controls[control_file] = controls
            return controls

class ControlVariants:
    def __init__(self, name, chrom, start, end, ref, control):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.ref = ref
        self.control = control
        self.exp_freq = float(name.split("_")[-1])
        self.variants = self._calc_diffs(ref, control)
        self._matches = set([])

    def get_no_matches(self):
        out = []
        matched = 0
        for vi, (_, ref, alt) in self.variants.items():
            if vi not in self._matches:
                out.append((vi + self.start, ref, alt))
            else:
                matched += 1
        print(self.name, self._matches, matched)
        return out

    def _calc_diffs(self, ref, control):
        """Calculate all differences between pre-aligned reference and control for validation.
        """
        vs = collections.OrderedDict()
        last_r, last_c = None, None
        cur_indel_r, cur_indel_c, cur_indel_i, cur_indel_ref_i = None, None, None, None
        ref_i = -1
        for i, r in enumerate(ref):
            if r != "-":
                ref_i += 1
            c = control[i]
            if c != r:
                if c == "-" or r == "-":
                    if not cur_indel_i:
                        cur_indel_r = last_r
                        cur_indel_c = last_c
                        cur_indel_i = i - 1
                        cur_indel_ref_i = ref_i - 1
                    cur_indel_r += r
                    cur_indel_c += c
                else:
                    if cur_indel_i:
                        vs[cur_indel_ref_i] = (cur_indel_i, cur_indel_r, cur_indel_c)
                        cur_indel_r, cur_indel_c, cur_indel_i, cur_indel_ref_i = None, None, None, None
                    vs[ref_i] = (i, r, c)
            elif cur_indel_i:
                vs[cur_indel_ref_i] = (cur_indel_i, cur_indel_r, cur_indel_c)
                cur_indel_r, cur_indel_c, cur_indel_i, cur_indel_ref_i = None, None, None, None
            last_r, last_c = r, c
        return vs

    def variant_matches(self, chrom, pos, ref_allele, allele):
        assert chrom == self.chrom
        lpos_ref = pos - self.start
        lpos = lpos_ref + self.ref[:lpos_ref].count("-")
        assert lpos >= 0

        cref = self.ref[lpos:lpos + len(ref_allele)]
        cref = self.ref[lpos:lpos + len(ref_allele) + cref.count("-")]
        calt = cref[:]
        alt_gaps = 0
        matches = []
        for ref_vi, (vi, vref, valt) in self.variants.items():
            if ref_vi >= lpos_ref and ref_vi < lpos_ref + len(ref_allele):
                vs = vi - lpos
                ve = vs + len(vref)
                if len(cref[vs:ve].replace("-", "")) == len(vref.replace("-", "")):
                    assert cref[vs:ve].replace("-", "") == vref.replace("-", ""), \
                        (cref[vs:ve], vref, len(cref), vi, vs, ve)
                    calt = calt[:vs - alt_gaps] + valt.replace("-", "") + calt[ve - alt_gaps:]
                    alt_gaps += valt.count("-")
                    matches.append(ref_vi)
        if calt == allele:
            for i in matches:
                self._matches.add(i)
            return True
        else:
            return False

def _prep_control_alignments(control_file, control_coords, ref_file):
    """Align controls to references to provide baseline for assessing variability.
    """
    ref = SeqIO.parse(ref_file, "fasta").next()
    chrom, start, end = control_coords
    assert chrom == ref.id
    ref_seq = str(ref.seq[start:end])
    aligns = []
    for rec in SeqIO.parse(control_file, "fasta"):
        rec_start, rec_end = [int(x) for x in rec.id.split("_")[-2].split("-")]
        start_offset = start - rec_start + 1
        end_offset = rec_end - end
        rec_seq = str(rec.seq)[start_offset:-end_offset]
        cur_align = pairwise2.align.globalms(ref_seq, rec_seq, 1, 0, -1, -0.1,
                                             one_alignment_only=True)[0][:2]
        aligns.append(ControlVariants(rec.id, chrom, start, end, cur_align[0], cur_align[1]))
    return aligns

def _check_haplotype_matches(cmanage, chrom, start, ref_allele, allele_freqs, ref_file,
                             region_file, config):

    for allele, freq, count in allele_freqs:
        if count > config["params"]["freebayes"]["min_count"]:
            control_matches = []
            for control in cmanage.get_controls(chrom, start):
                if control.variant_matches(chrom, start, ref_allele, allele):
                    control_matches.append(control)
            if len(control_matches) == 0:
                print("fp", chrom, start, allele, freq, count)
                cmanage.evals["fps"]["%.1f" % freq] += 1
            else:
                exp_freq = sum(m.exp_freq for m in control_matches)
                if abs(1 - (freq / exp_freq)) < 0.4 or abs(freq - exp_freq) < 1.0:
                    cmanage.evals["tps"]["%.1f" % exp_freq] += 1
                else:
                    print("wrongfreq", chrom, start, allele, freq, exp_freq)
                    cmanage.evals["wrongfreq"]["%.1f" % exp_freq] += 1
    return cmanage

# -- Variant calling

def _run_vardict(bam_file, ref_file, region_file, config):
    """Variant calling with VarDict
    """
    work_dir = "vardict"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    out_file = os.path.join(work_dir, "%s.vcf" % os.path.splitext(os.path.basename(bam_file))[0])
    if not os.path.exists(out_file):
        with pysam.Samfile(bam_file) as in_bam:
            name = in_bam.header["RG"][0]["SM"]
        min_af = config["params"]["vardict"]["min_af"]
        cmd = ("vardict-java -G {ref_file} -f {min_af} -N {name} -b {bam_file} "
            "-c 1 -S 2 -E 3 -g 4 {region_file} | "
            "teststrandbias.R | var2vcf_valid.pl -N {name} -f {min_af} | "
            "vcfallelicprimitives | vt normalize -r {ref_file} - > {out_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    filter_file = "%s-filter%s" % os.path.splitext(out_file)
    with open(out_file) as in_handle:
        with open(filter_file, "w") as out_handle:
            for line in in_handle:
                remove = False
                if not line.startswith("#"):
                    if "pSTD" in line.split("\t")[6].split(";"):
                        remove = True
                if not remove:
                    out_handle.write(line)
    return filter_file

def _freebayes_call(bam_file, ref_file, region_file, config):
    """Run FreeBayes calling, trying to emphasize
    """
    out_file = os.path.join("freebayes", "%s.vcf" % os.path.splitext(os.path.basename(bam_file))[0])
    if not os.path.exists(out_file):
        min_alt_fraction = config["params"]["freebayes"]["min_alt_fraction"]
        hl = config["params"]["freebayes"]["haplotype_length"]
        haplotype_length = "--haplotype-length %s" % (hl) if hl > 0 else ""
        cmd = ("freebayes {bam_file} -f {ref_file} --pooled-continuous "
               "--pooled-discrete --report-genotype-likelihood-max --allele-balance-priors-off "
               "--min-alternate-fraction {min_alt_fraction} --min-alternate-count 2 "
               "{haplotype_length} -t {region_file} | "
               "vcfallelicprimitives | vt normalize -r {ref_file} - | vcfuniqalleles "
               "> {out_file}")
        print(cmd.format(**locals()))
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

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

def _check_pair_overlap(fq1_file, fq2_file, ref_file, config):
    """Ensure that pairs of fastq files overlap
    """
    import matplotlib.pyplot as plt
    check_file = "%s-check_ol.bam" % (fq1_file)
    if not os.path.exists(check_file):
        cmd = ("bwa mem -M {ref_file} {fq1_file} {fq2_file} | "
               "samtools sort - -n -o {check_file} -O bam -T {check_file}-tmp")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    overlaps = 0
    total = 0
    starts = collections.defaultdict(int)
    ends = collections.defaultdict(int)
    with pysam.Samfile(check_file) as in_bam:
        for _, pair in itertools.groupby(in_bam, key=lambda x: x.query_name):
            coords = []
            for r in pair:
                if not r.is_secondary and not r.is_unmapped:
                    coords.append((r.reference_start, r.reference_end))
            total += 1
            if len(coords) == 2:
                if len(set(range(*coords[0])) & set(range(*coords[1]))) > config["params"]["pear"]["min_overlap"]:
                    overlaps += 1
                    starts[min([x[0] for x in coords])] += 1
                    ends[max([x[1] for x in coords])] += 1
    print(overlaps, total, "%.1f%%" % (float(overlaps) / total * 100.0))
    starts = sorted(list(starts.items()))
    ends = sorted(list(ends.items()))
    plt.clf()
    plt.plot([x[0] for x in starts], [x[1] for x in starts], [x[0] for x in ends], [x[1] for x in ends])
    plt.savefig("%s.png" % os.path.splitext(check_file)[0])

def _plot_coverage(bam_file):
    """Simple coverage plot for an input BAM file.
    """
    import matplotlib.pyplot as plt
    out_file = "%s-coverage.pdf" % os.path.splitext(bam_file)[0]
    if not os.path.exists(out_file):
        if not os.path.exists(bam_file + ".bai"):
            subprocess.check_call(["samtools", "index", bam_file])
        x, y = [], []
        with pysam.Samfile(bam_file) as in_bam:
            for pileup in in_bam.pileup(max_depth=1e6):
                if pileup.nsegments > 10:
                    x.append(pileup.reference_pos)
                    y.append(pileup.nsegments)
        plt.clf()
        plt.plot(x, y)
        plt.savefig(out_file)
    return out_file

def _realign_merged(fq1_file, fq2_file, region_file, ref_file):
    """Realign merged reads back to reference genome.
    """
    out_file = "%s.bam" % os.path.splitext(fq1_file)[0]
    if not fq2_file:
        fq2_file = ""
    if region_file:
        region_cmd = " | bedtools intersect -abam /dev/stdin -b %s" % region_file
    else:
        region_cmd = ""
    if not os.path.exists(out_file):
        cmd = ("bwa mem {ref_file} {fq1_file} {fq2_file} | "
               "samtools sort - -O bam -T {out_file}-tmp " +
               region_cmd +
               " > {out_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    if not os.path.exists(out_file + ".bai"):
        subprocess.check_call(["samtools", "index", out_file])
    return out_file

def _select_full_fastqs(bam_file, work_dir):
    base, ext = os.path.splitext(os.path.basename(bam_file))
    base = os.path.join(work_dir, base)
    full_fq1 = "%s-1.fastq" % base
    full_fq2 = "%s-2.fastq" % base
    if not os.path.exists(full_fq1):
        cmd = ("samtools view -u -f 0x2 {bam_file} | samtools sort -n -O bam -T {full_fq1} - | "
               "bedtools bamtofastq -i /dev/stdin -fq {full_fq1} -fq2 {full_fq2}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return full_fq1, full_fq2

def _select_fastq_in_region(bam_file, region, region_file):
    """Extract paired end reads where either read overlaps the region of interest.
    """
    full_fq1, full_fq2 = _select_full_fastqs(bam_file)
    base, ext = os.path.splitext(bam_file)

    chrom, start, end = region.split("-")
    fq1 = "%s-%s-1.fastq" % (base, region)
    fq2 = "%s-%s-2.fastq" % (base, region)
    name_file = "%s-%s-names.txt" % (base, region)
    keep_file = "%s-keep%s" % os.path.splitext(name_file)
    if not os.path.exists(keep_file):
        cmd = ("bedtools intersect -abam {bam_file} -b {region_file} | samtools view - "
               "| cut -f 1 | sort | uniq > {name_file}")
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
    return [fq1, fq2]

def _subset_region_file(in_file):
    """Subset a BED file into individual regions to be analyzed separately.
    """
    out_dir = os.path.join(os.getcwd(), "regions")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out = []
    with open(in_file) as in_handle:
        for line in in_handle:
            chrom, start, end, name = line.strip().split("\t")[:4]
            base, ext = os.path.splitext(os.path.basename(in_file))
            out_file = os.path.join(out_dir, "%s-%s-%s-%s%s" % (base, chrom, start, end, ext))
            with open(out_file, "w") as out_handle:
                out_handle.write(line)
            out.append(("%s-%s-%s" % (chrom, start, end), name, out_file))
    return out

# reference management

def _subset_viquas_file(in_file, region, ref_file):
    """Subset output file by regions that differ from the reference genome.

    Identify those matching to reference by ends and write regions of
    difference from reference.
    """
    pad = 15
    end_buffer = 100
    base, ext = os.path.splitext(in_file)
    out_file = "%s-%s%s" % (base, region, ext)
    out_file_bad = "%s-nomatch%s" % os.path.splitext(out_file)
    chrom, start, end = region.split("-")
    orig_start = int(start) - 1
    orig_end = int(end)
    start = orig_start - end_buffer
    end = orig_end + end_buffer
    ref_index = pyfaidx.Fasta(ref_file)
    ref_seq = ref_index[chrom][start:end].seq
    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        index = pyfaidx.Fasta(in_file)
        with open(out_file, "w") as out_handle:
            with open(out_file_bad, "w") as out_handle_bad:
                for recid in index.keys():
                    cur_seq = index[recid][:].seq
                    cur_start = cur_seq.find(ref_seq[:pad])
                    added = False
                    if cur_start > 0:
                        cur_end = cur_seq.find(ref_seq[-pad:], cur_start)
                        if cur_end > 0:
                            added = True
                            cur_seq = cur_seq[cur_start:cur_end + pad]
                            cur_align = pairwise2.align.globalms(cur_seq, ref_seq, 1, 0, -1, -0.1,
                                                                 one_alignment_only=True)
                            match_start, match_end = _find_match_aligns(cur_align[0][0], cur_align[0][1], pad,
                                                                        end_buffer)
                            final_start = match_start + cur_start
                            final_end = match_end + cur_start
                            rec = index[recid][final_start:final_end]
                            out_handle.write(repr(rec) + "\n")
                    if not added:
                        out_handle_bad.write(repr(index[recid][orig_start:orig_end]) + "\n")
    return out_file

def _find_match_aligns(seq1, seq2, pad, end_buffer):
    min_ends = end_buffer + pad * 2
    for s_i, char2 in enumerate(seq2):
        if s_i > min_ends:
            if seq1[s_i] != char2:
                break
    seq1_r = list(reversed(list(seq1)))
    for e_i, char2 in enumerate(reversed(list(seq2))):
        if e_i > min_ends:
            if seq1_r[e_i] != char2:
                break
    return s_i - pad, len(seq1.replace("-", "")) - e_i + pad

def _subset_ref_file(ref_file, out_dir, region):
    chrom, start, end = region.split("-")
    base, ext = os.path.splitext(os.path.basename(ref_file))
    out_file = os.path.join(out_dir, "%s-%s%s" % (base, region, ext))
    region = "%s:%s-%s" % (chrom, int(start), int(end))
    cmd = "samtools faidx {ref_file} {region} | sed 's/{region}/{chrom}/' > {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def _index_ref_file(ref_file):
    out_file = "%s.fai" % ref_file
    if not os.path.exists(out_file):
        cmd = ["samtools", "faidx", ref_file]
        subprocess.check_call(cmd)
    bwa_index = "%s.bwt" % ref_file
    if not os.path.exists(bwa_index):
        cmd = ["bwa", "index", ref_file]
        subprocess.check_call(cmd)
    return out_file

# control library evaluation

def _is_control(hap_file):
    base = os.path.basename(os.path.dirname(hap_file))
    return base.startswith(("1_", "2_")) or "Control" in base

def _evaluate_control(hap_file, ref_file):
    thresh = 0.9
    thresh_span = 0.7
    hap_file = os.path.relpath(hap_file, os.getcwd())
    cmd = ["ssearch36", hap_file, ref_file, "-d", "1", "-m", "9i", "-m", "4"]
    output = subprocess.check_output(cmd)
    cur_query = None
    query_length = None
    in_best = False
    results = []
    for line in output.split("\n"):
        if line.strip().startswith("%s>>>" % (len(results) + 1)):
            parts = line.split(" - ")
            cur_query = parts[0].strip().split(">>>")[1]
            if len(parts) > 1:
                query_length = float(parts[1].split()[-2].strip())
        elif cur_query and line.startswith("The best scores are:"):
            in_best = True
        elif cur_query and in_best:
            if not line:
                print("Did not find match for", cur_query)
                in_best = False
                cur_query = None
            else:
                info = line.split()
                hit = info[0]
                match_pct = float(info[-3])
                span_pct = float(info[-1]) / query_length
                if span_pct > 1.1:
                    print(query_length, line)
                    raise NotImplementedError
                if span_pct > thresh_span:
                    if match_pct > thresh:
                        results.append((cur_query, hit, match_pct, span_pct))
                        in_best = False
                        cur_query = None
    ctl_sum = collections.defaultdict(float)
    eval_file = "%s-eval.txt" % os.path.splitext(hap_file)[0]
    with open(eval_file, "w") as eval_handle:
        for orig, match, match_pct, span_pct in results:
            print(orig, match, match_pct, span_pct, file=eval_handle)
            _, freq = orig.split(":")[0].split("_")
            ctl_sum[match] += float(freq)
        pprint.pprint(dict(ctl_sum), stream=eval_handle)

def compare_vcf(config_file, vcf_file, region_file):
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    _summarize_calls(vcf_file, config["ref_file"], config["control_vcf"],
                     region_file, config["params"]["validation"])

if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[2].endswith(".vcf") and sys.argv[3].endswith(".bed"):
        compare_vcf(*sys.argv[1:])
    else:
        main(*sys.argv[1:])
