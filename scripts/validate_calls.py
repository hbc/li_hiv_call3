#!/usr/bin/env python
"""Validate VCF file of calls against a truth VCF.

Usage:
  validate_calls.py calls.vcf control.vcf regions.bed params.yaml
"""
import collections
import sys

import numpy
import pysam
import yaml

def main(call_vcf, control_vcf, region_bed, param_file):
    with open(param_file) as in_handle:
        params = yaml.safe_load(in_handle)
    _summarize_calls(call_vcf, control_vcf, region_bed, params["params"]["validation"])

def _summarize_calls(vcf_file, control_file, region_file, params):
    """Summarize calls against standard validation VCF, not requiring phasing.
    """
    want_pos = _read_regions(region_file, params["end_buffer"])
    control_calls = _read_control(control_file, want_pos)

    tp_reported = collections.defaultdict(list)
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
                        tp_reported[control_calls[key]].append(freq)
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
    print("** TP/FP/FN by frequency")
    freqs = list(sorted(set(tp_freqs.keys() + fn_freqs.keys()), reverse=True))
    fp_freqs = collections.defaultdict(int)
    for fp in sorted(fps_lowfreq + fps):
        cur_freq = float(fp[1]["AF"]) * 100.0
        freq = min(freqs, key=lambda x: abs(x - cur_freq))
        fp_freqs[freq] += 1
    print("| freq | TP | FN | FP |")
    for freq in freqs:
        print("| %.1f | %s | %s | %s |" % (freq, tp_freqs[freq], fn_freqs[freq], fp_freqs[freq]))
    print "** Called frequencies for TPs"
    print("| freq | median | min | max |")
    for freq in freqs:
        called_freqs = tp_reported[freq]
        if len(called_freqs):
            print("| %.1f | %.1f | %.1f | %.1f |" % (freq, numpy.median(called_freqs),
                                                     min(called_freqs), max(called_freqs)))
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

if __name__ == "__main__":
    main(*sys.argv[1:])
