import sys

import os
from os import listdir
from os.path import isfile, isdir, join


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import gffutils

import multiprocessing

import edlib
import sys
import re

def load_fasta(filename, tp = "list"):
    if tp == "map":
        record_lst = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    else:
        record_mp = list(SeqIO.parse(filename, "fasta"))
        record_lst = {}
        for a in record_mp:
            record_lst[a.id] = a
    return record_lst

def edist(lst):
    if len(str(lst[0])) == 0:
        return len(str(lst[1]))
    if len(str(lst[1])) == 0:
        return len(str(lst[0]))
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    return result["editDistance"], result["cigar"]

def aai(ar):
    p1, p2 = ar[0], ar[1]
    if p1.endswith("*"):
        p1 = p1[:-1]
    ed, cigar = edist([str(p1), str(p2)])
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= max(len(p1), len(p2))
    return aai

def load_tsv(filename, uniprot, synth):
    res = set()
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            qseqid, sseqid, pident, qlen, slen, length, qstart, qend, sstart, send = ln.strip().split("\t")
            if qseqid not in res and float(pident) >= 90 and int(length) > 0.8*int(slen) and int(length) > 0.8*int(qlen):
                if aai([uniprot[qseqid].seq, synth[sseqid].seq]) > 0.9:
                    res.add(qseqid)
    return res

alignments = load_fasta(sys.argv[1], "map")
uniprot = load_fasta(sys.argv[2], "map")
synth = load_fasta(sys.argv[3], "map")
mapped = load_tsv(sys.argv[2], uniprot, synth)

total_mapped = 0
has_alignment_tp = 0
has_alignment_fp = 0
has_alignment_tn =0
has_alignment_fn =0
for u in uniprot:
    if u in mapped:
        total_mapped += 1
        if "SPAligner_" + u in alignments:
            has_alignment_tp += 1
        else:
            has_alignment_fn += 1
    else:
        if "SPAligner_" + u in alignments:
            has_alignment_fp += 1
        else:
            has_alignment_tn += 1

print "total_mapped", total_mapped
print "tp", has_alignment_tp, "fp", has_alignment_fp, "tn", has_alignment_tn, "fn", has_alignment_fn





