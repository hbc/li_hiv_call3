"""Create VCF of control population from input FASTA.

Explore approaches to cleanly create a single truth VCF with multiple
frequency information.
"""
import os
import subprocess
import sys

from Bio import SeqIO

def main(control_file, ref_file, jvarkit_path=None):
    approach = "pileup"
    start, end = _get_coords(control_file)
    work_dir = os.path.join(os.getcwd(), "%s-prep" % os.path.splitext(os.path.basename(control_file))[0])
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if approach == "pileup":
        pileup_call_all(control_file, ref_file, work_dir)
    elif approach == "multialign":
        contig, input_file = make_input_file(control_file, start, end, ref_file, work_dir)
        out_vcf = prep_vcf(input_file, contig, start, jvarkit_path, ref_file)
        print out_vcf

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

def pileup_call_all(control_file, ref_file, work_dir):
    """Call pileup-based variant calling for a set of inputs in a control file.
    """
    for rec in SeqIO.parse(control_file, "fasta"):
        out_file = os.path.join(work_dir, "%s.fa" % rec.description.replace(".", "_"))
        with open(out_file, "w") as out_handle:
            SeqIO.write([rec], out_handle, "fasta")
        call_with_alignment(out_file, ref_file)

def call_with_alignment(in_file, ref_file):
    print in_file

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
