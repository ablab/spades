import sys
import pysam
import subprocess
import matplotlib
import edlib

from Bio import SeqIO

reads_fasta = sys.argv[1]
reference_fasta = sys.argv[2]
grap_path = sys.argv[3]
K = sys.argv[4]
tp = sys.argv[5]
idealreads_aligner = sys.argv[6]

if len(sys.argv) < 7:
    print "Builds ideal paths.\nUsage: python2 build_idealpaths.py reads.fasta reference.fa graph_saves reads_type{pacbio, nanopore} ideal_seq_aligner_path"
    exit(-1)

def edist(lst):
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", additionalEqualities=[('U', 'T')
                                                , ('R', 'A'), ('R', 'G')
                                                , ('Y', 'C'), ('Y', 'T'), ('Y', 'U')
                                                , ('K', 'G'), ('K', 'T'), ('K', 'U')
                                                , ('M', 'A'), ('M', 'C')
                                                , ('S', 'C'), ('S', 'G')
                                                , ('W', 'A'), ('W', 'T'), ('W', 'U')
                                                , ('B', 'C'), ('B', 'G'), ('B', 'T'), ('B', 'U')
                                                , ('D', 'A'), ('D', 'G'), ('D', 'T'), ('D', 'U')
                                                , ('H', 'A'), ('H', 'C'), ('H', 'T'), ('H', 'U')
                                                , ('V', 'A'), ('V', 'C'), ('V', 'G')
                                                , ('N', 'C'), ('N', 'C'), ('N', 'G'), ('N', 'T'), ('N', 'U')] )
    return result["editDistance"]   

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
    mp_real = {}
    mp_edlib = {}
    cnt = 0
    ed = {}
    for r in readssam:
        start = r.reference_start
        end = r.reference_end
        name = r.query_name
        if start == -1 or end is None:
            continue
        sequence = reffasta.fetch(r.reference_name, start, end)
        aligned_sequence = r.query_alignment_sequence
        cnt += 1
        if (name not in mp or len(mp[name]) < len(sequence)) and len(sequence) > 0.8*len(reads[name].seq):
            if r.is_reverse:
                sequence = make_rc(sequence)
                aligned_sequence = make_rc(aligned_sequence)
            mp[name] = sequence
            mp_real[name] = aligned_sequence
            ed[name] = edist([sequence, aligned_sequence])*1.0/len(sequence) 
    readssam.close()
    reffasta.close()
    avg_identity = 0.0
    for it in ed.keys():
        avg_identity += ed[it]
    print("Number of mapped items > 80%: " + str(len(mp)))
    print("Number of mapped items: " + str(cnt))
    print("Number of items: " + str(len(reads)))
    print("Avg identity (> 80%): " + str(1 - avg_identity/len(ed)))
    return [mp, mp_real]


def save(refseq_fasta, seq):
    with open(refseq_fasta, "w") as out:
        for k in seq.keys():
            out.write(">" + k + "\n" + seq[k] + "\n")


def run_bwamem(reads_fasta, reference_fasta, tp, build_index = False):
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
    [refseq, realseq] = load_ideal_reads(reads_bam, reference_fasta, reads)
    refseq_fasta =  "/".join(reads_fasta.split("/")[:-1] + [ "refseq_" + reads_fasta.split("/")[-1]])
    print("> Save into " + refseq_fasta)
    save(refseq_fasta, refseq)

    realseq_fasta =  "/".join(reads_fasta.split("/")[:-1] + [ "realseq_" + reads_fasta.split("/")[-1]])
    print("> Save into " + realseq_fasta)
    save(realseq_fasta, realseq)

    return refseq_fasta

def run_minimap2(reads_fasta, reference_fasta, tp, build_index = False):
    print("> Load reads from fasta")
    reads = load_fasta(reads_fasta)
    reads_bam = reads_fasta[:-len(".fasta")] + "_origin.bam"
    print("> Run minimap2 to align reads on reference and save result to " + reads_bam)
    if tp == "pacbio":
       t = "map-pb"
    elif tp == "nanopore":
       t = "map-ont"
    print("minimap -ax " + t + " -t 33 " + reference_fasta + " " + reads_fasta + " | samtools view -bS - > " + reads_bam)
    subprocess.call(["minimap -ax " + t + " -t 33 " + reference_fasta + " " + reads_fasta + " | samtools view -bS - > " + reads_bam], shell=True)
    print("> Load mappings from bam (longest mapping for each read name and not less 20% length)")
    refseq = load_ideal_reads(reads_bam, reference_fasta, reads)
    refseq_fasta =  "/".join(reads_fasta.split("/")[:-1] + [ "refseq_minimap" + reads_fasta.split("/")[-1]])
    print("> Save into " + refseq_fasta)
    save(refseq_fasta, refseq)
    return refseq_fasta

def run_idealreads_aligner(refseq_fasta):
    print("> Run reference sequence aligner ")
    print(idealreads_aligner + " " + K + " " + grap_path + " " + refseq_fasta + " " + refseq_fasta[:-6] + "_mapping --spades > " + refseq_fasta[:-6] + ".log")
    subprocess.call([idealreads_aligner + " " + K + " " + grap_path + " " + refseq_fasta + " " + refseq_fasta[:-6] + "_mapping --spades > " + refseq_fasta[:-6] + ".log"], shell=True)


refseq_fasta = run_bwamem(reads_fasta, reference_fasta, tp)
run_idealreads_aligner(refseq_fasta)


