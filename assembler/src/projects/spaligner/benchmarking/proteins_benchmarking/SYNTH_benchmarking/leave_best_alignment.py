from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

import edlib 
import sys
import re

def edist(lst):
    if len(str(lst[0])) == 0:
        return len(str(lst[1]))
    if len(str(lst[1])) == 0:
        return len(str(lst[0]))
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    return result["editDistance"], result["cigar"]

def aai(p1, p2):
    if p1.endswith("*"):
        p1 = p1[:-1]
    ed, cigar = edist([str(p1), str(p2)])
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= max(len(p1), len(p2))
    return aai

def load_fasta(filename, tp = "list"):
    if tp == "map":
        record_lst = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    else:
        record_lst = list(SeqIO.parse(filename, "fasta"))
    return record_lst 


alignments = load_fasta(sys.argv[1])
proteins = load_fasta(sys.argv[2], "map")
best_alignments = {}
best_score = {}
with open(sys.argv[1][:-len(".fasta")] + "_best.fasta", "w") as fout:
    for n in alignments:
        name = n.name
        lst = n.description.split("|")
        start_s, end_s, score = int(lst[-3][len("start_s="):]), int(lst[-2][len("end_s="):]), int(lst[-1][len("score="):])
        aai_v = aai(n.seq.translate(), proteins[name].seq) 
        if aai_v >= 0.9 and (name not in best_score or best_score[name] < aai_v):
            best_alignments[name] = str(n.seq.translate())
            best_score[name] = aai_v
    for al in best_alignments:
        fout.write(">SPAligner_" + al + "\n" + best_alignments[al] + "\n")
