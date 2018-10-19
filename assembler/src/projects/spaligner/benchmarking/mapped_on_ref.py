import sys
import pysam
import subprocess
import matplotlib

from Bio import SeqIO

reads_fasta = sys.argv[1]
tp = sys.argv[2]
reference_fasta = sys.argv[3]

def load_fasta(filename):
    record_dict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    return record_dict

def make_rc(seq):
    mapping = {"A": "T", "T": "A", "C": "G", "G": "C"}
    res = []
    for c in seq:
        res.append(mapping[c])
    return "".join(res[::-1])

def load_ideal_reads(reads_bam, reference_fasta, reads):
    readssam = pysam.AlignmentFile(reads_bam, "rb")
    reffasta = pysam.FastaFile(reference_fasta)

    mp = {}
    for r in readssam:
        start = r.reference_start
        end = r.reference_end
        name = r.query_name
        if start == -1 or end is None:
            continue
        sequence = reffasta.fetch(r.reference_name, start, end)
        if sequence == None:
            continue
        if (name not in mp or len(mp[name]) < len(sequence)) and len(sequence) > 0.8*len(reads[name].seq):
            if r.is_reverse:
                sequence = make_rc(sequence)
            mp[name] = sequence
    readssam.close()
    reffasta.close()
    print("Number of mapped items > 80%: " + str(len(mp)))
    print("Number of items: " + str(len(reads)))

def run_minimap2(reads_fasta, reference_fasta, tp):
    print("> Load reads from fasta")
    reads = load_fasta(reads_fasta)
    reads_bam = reads_fasta[:-len(".fasta")] + "_origin.bam"
    print("> Run minimap2 to align reads on reference and save result to " + reads_bam)
    if tp == "pacbio":
       t = "map-pb"
    elif tp == "nanopore":
       t = "map-ont"
    print("/home/tdvorkina/soft/minimap2/minimap2/minimap2 -ax " + t + " -t 33 " + reference_fasta + " " + reads_fasta + " | samtools view -bS - > " + reads_bam)
    subprocess.call(["/home/tdvorkina/soft/minimap2/minimap2/minimap2 -ax " + t + " -t 33 " + reference_fasta + " " + reads_fasta + " | samtools view -bS - > " + reads_bam], shell=True)
    print("> Load mappings from bam (longest mapping for each read name and not less 20% length)")
    load_ideal_reads(reads_bam, reference_fasta, reads)


def run_bwamem(reads_fasta, reference_fasta, tp):
    print("> Load reads from fasta")
    reads = load_fasta(reads_fasta)
    reads_bam = reads_fasta[:-len(".fasta")] + "_origin.bam"
    print("> Run bwa mem to align reads on reference and save result to " + reads_bam)
    if tp == "pacbio":
        t = "pacbio"
    elif tp == "nanopore":
        t = "ont2d"
    print("bwa mem -x " + t +  " -t 33 " + reference_fasta + " " + reads_fasta + " | samtools view -bS - > " + reads_bam)
    subprocess.call(["bwa mem -x " + t + " -t 33 " + reference_fasta + " " + reads_fasta + " | samtools view -bS - > " + reads_bam], shell=True)
    print("> Load mappings from bam (longest mapping for each read name and not less 20% length)")
    load_ideal_reads(reads_bam, reference_fasta, reads)

run_minimap2(reads_fasta, reference_fasta, tp)
run_bwamem(reads_fasta, reference_fasta, tp)