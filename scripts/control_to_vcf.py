"""Create VCF of control population from input FASTA.

Explore approaches to cleanly create a single truth VCF with multiple
frequency information.
"""
import collections
import copy
import os
import subprocess
import sys

from Bio import SeqIO
import toolz as tz
import vcf
from vcf import model as vcf_model
from vcf import parser as vcf_parser
import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    if tz.get_in(["params", "validation", "vcf_prep_method"], config) == "pileup":
        out_vcf = pileup_prep_controls(config["controls"], config["ref_file"])
    elif tz.get_in(["params", "validation", "vcf_prep_method"], config) == "multialign":
        contig, input_file = make_input_file(control_file, start, end, ref_file, work_dir)
        out_vcf = prep_vcf(input_file, contig, start, jvarkit_path, ref_file)
    else:
        out_vcf = None
    print "Control VCF", out_vcf

def _get_coords(in_file):
    """Retrieve coordinates from the headers of input FASTA files.
    """
    coords = []
    with open(in_file) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                coord_str = line.split("_")[-2]
                start, end = coord_str.split("-")
                coords.append(int(start))
                coords.append(int(end))
    return min(coords), max(coords)

# -- pileup based variant calling

def pileup_prep_controls(control_info, ref_file):
    vcfs = []
    for name, control_file in control_info.items():
        start, end = _get_coords(control_file)
        work_dir = os.path.join(os.getcwd(), "%s-prep" % os.path.splitext(os.path.basename(control_file))[0])
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
        vcf_file = pileup_call_control(control_file, ref_file, work_dir)
        vcfs.append(((start, end), vcf_file))
    vcfs = [xs[-1] for xs in sorted(vcfs)]
    combined_calls = concat_control_vcfs(vcfs, os.getcwd())
    return calculate_freqs(combined_calls)

def calculate_freqs(in_file):
    CallData = collections.namedtuple('calldata', ["GT"])
    out_file = "%s-freqs.vcf" % in_file.replace(".vcf.gz", "")
    in_bcf = vcf.Reader(filename=in_file)
    template = copy.copy(in_bcf)
    template.samples = ["control"]
    template.infos["FREQ"] = vcf_parser._Info("FREQ", "A", "Float", "Control frequencies for each allele", "", "")
    out_bcf = vcf.Writer(open(out_file, "w"), template)
    print out_bcf.template.infos
    for rec in in_bcf:
        freqs = {}
        for i in range(len(rec.ALT)):
            freqs[i + 1] = 0.0
        for sample in rec.samples:
            for ai in sample.data.GT.split("/") if sample.data.GT.find("/") > 0 else sample.data.GT.split("|"):
                if ai != "." and int(ai) > 0:
                    freqs[int(ai)] += float(sample.sample)
        freqs = ",".join(["%.1f" % freq for i, freq in sorted(freqs.items())])
        rec.INFO = {}
        rec.add_info("FREQ", freqs)
        rec.samples = [vcf_model._Call(rec, "control", CallData("0/1"))]
        out_bcf.write_record(rec)
    return out_file

def concat_control_vcfs(vcfs, work_dir):
    out_file = os.path.join(os.getcwd(), "controls.vcf.gz")
    in_vcf_str = " ".join(vcfs)
    cmd = "bcftools concat -O z -o {out_file} {in_vcf_str}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    cmd = ("tabix -f -p vcf {out_file}")
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def pileup_call_control(control_file, ref_file, work_dir):
    """Call pileup-based variant calling for a set of inputs in a control file.
    """
    inputs = []
    for rec in SeqIO.parse(control_file, "fasta"):
        out_file = os.path.join(work_dir, "%s.fa" % rec.description.replace(".", "_"))
        rec.description = rec.description.split("_")[-1]
        with open(out_file, "w") as out_handle:
            SeqIO.write([rec], out_handle, "fasta")
        inputs.append((float(rec.description), out_file))
    inputs.sort(reverse=True)
    vcfs = []
    for freq, fasta_file in inputs:
        vcfs.append(call_with_alignment(fasta_file, ref_file, "%.1f" % freq))
    return merge_pileup_calls(vcfs, control_file, ref_file, work_dir)

def call_with_alignment(in_file, ref_file, sample_name):
    """Do alignment and pileup based calling to produce a VCF file for a sample.
    """
    align_file = "%s.bam" % os.path.splitext(in_file)[0]
    cmd = (r"bwa mem -R '@RG\tID:{sample_name}\tPL:illumina\tPU:{sample_name}\tSM:{sample_name}' "
           r"{ref_file} {in_file}"
           r" | samtools sort -O bam -T {align_file}-tmp -o {align_file}")
    subprocess.check_call(cmd.format(**locals()), shell=True)
    variant_file = "%s.vcf.gz" % os.path.splitext(in_file)[0]
    cmd = ("freebayes -f {ref_file} {align_file} --min-alternate-fraction 0 --pooled-continuous "
           "--report-monomorphic --haplotype-length 0 --min-alternate-count 1 | "
           "vcfallelicprimitives | vt normalize -r {ref_file} - | "
           "sed 's:0/0:0/1:' | sed 's:0|0:0|1:' | "
           "bgzip -c > {variant_file}")
    subprocess.check_call(cmd.format(**locals()), shell=True)
    cmd = ("tabix -f -p vcf {variant_file}")
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return variant_file

def merge_pileup_calls(vcfs, control_file, ref_file, work_dir):
    out_file = os.path.join(work_dir, "%s.vcf.gz" % os.path.splitext(os.path.basename(control_file))[0])
    input_vcf_str = " ".join(vcfs)
    cmd = "bcftools merge -o {out_file} -O z {input_vcf_str}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

# -- multiple alignment based

def prep_vcf(input_file, contig, start, jvarkit_path, ref_file):
    aln_file = "%s-out.fa" % os.path.splitext(input_file)[0]
    cmd = "clustalo -i {input_file} -o {aln_file} --force"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    vcf_file = "%s.vcf" % os.path.splitext(aln_file)[0]
    cmd = "java -jar {jvarkit_path}/dist-*/biostar94573.jar -R {contig} {aln_file} > {vcf_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    fix_file = "%s-fix%s" % os.path.splitext(vcf_file)
    with open(vcf_file) as in_handle:
        with open(fix_file, "w") as out_handle:
            for line in in_handle:
                if not line.startswith("#"):
                    parts = line.split("\t")
                    parts[1] = str(int(parts[1]) + start - 1)
                    line = "\t".join(parts)
                out_handle.write(line)
    norm_file = "%s-norm%s" % os.path.splitext(fix_file)
    cmd = "vt normalize -r {ref_file} {fix_file} | bcftools norm -f {ref_file} --check-ref w - > {norm_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return norm_file

def make_input_file(control_file, start, end, ref_file, work_dir):
    with open(ref_file) as in_handle:
        contig = in_handle.readline().strip().replace(">", "")
    out_file = os.path.join(work_dir, os.path.basename(control_file))
    cmd = "samtools faidx {ref_file} {contig}:{start}-{end} | sed 's/{contig}:{start}-{end}/{contig}/' > {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    cmd = "cat {control_file} >> {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return contig, out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
