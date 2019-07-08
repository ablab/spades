from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

import edlib 
import sys
import re


PATH2ARGS_FASTA =
PATH2SPALIGNER_RES =
PATH2CONTIGS =
PATH2AMR_RES =

def load_fasta(filename, tp = "list"):
    if tp == "map":
        record_lst = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    else:
        record_lst = list(SeqIO.parse(filename, "fasta"))
    return record_lst 

def edist(lst):
    if len(str(lst[0])) == 0:
        return len(str(lst[1]))
    if len(str(lst[1])) == 0:
        return len(str(lst[0]))
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    return result["editDistance"], result["cigar"]

def aai(ar):
    p1, p2 = str(ar[0]), str(ar[1])
    if p1.endswith("*"):
        p1 = p1[:-1]
    if p2.endswith("*"):
        p2 = p2[:-1]
    ed, cigar = edist([str(p1), str(p2)])
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= max(len(p1), len(p2))
    return aai*100 

def extract_amrfinder_alignments(filename, contigs, min_len, genes, prefix = "AMR_"):
    res = []
    with open(filename, "r") as fin:
        fin.readline()
        for ln in fin.readlines():
            lst = ln.strip().split("\t")
            name, start, end, strand, len_prop, acc_num = lst[0], int(lst[2]), int(lst[3]), lst[4], float(lst[10]), prefix + lst[13]
            if len_prop >= min_len:
                contig = contigs[name].seq
                if strand == "-":
                    if len(contig[start-1:end].reverse_complement()) % 3 != 0:
                        print name, start, end, strand, len_prop, acc_num
                    alignment = str(contig[start-1:end].reverse_complement().translate())
                else:
                    alignment = str(contig[start-1:end].translate())
            else:
                continue
            if aai([alignment, genes[lst[13]].seq]) >= identity:
                res.append({"name": prefix + lst[13], "aln": alignment})
    return res

def extract_spaligner_alignments(filename, min_len, identity, genes):
    alignments = load_fasta(filename)
    best_alignments = {}
    best_score = {}
    for n in alignments:
        name = n.name
        lst = n.description.split("|")
        start_s, end_s = int(lst[-4][len("start_s="):]), int(lst[-3][len("end_s="):])
        n_t = n.seq.translate()
        c_aai = aai([n_t, genes[name].seq])
        if len(str(n_t)) >= min_len/100.0*len(genes[name].seq) and c_aai >= identity and (name not in best_alignments or best_score[name] < c_aai):
            best_alignments[name] = str(n_t)
            best_score[name] = c_aai
    res = []
    for al in best_alignments:
        res.append({"name": "SP_" + al, "aln": best_alignments[al]})
    return res

assembler = "spades"
identity = 90
min_len = 0

proteins = load_fasta(PATH2ARGS_FASTA, "map")

input_data = {"spades": {"aln": PATH2SPALIGNER_RES, "contigs": PATH2CONTIGS, "amr": PATH2AMR_RES}}

contigs = load_fasta(input_data["spades"]["contigs"], "map")

alignments = extract_amrfinder_alignments(input_data[assembler]["amr"], contigs, 75, proteins)
print "Amr num", len(alignments)
was = len(alignments)
alignments.extend(extract_spaligner_alignments(input_data[assembler]["aln"], min_len, identity, proteins))
print "SPAligner num", len(alignments) - was



alignments_mp = {}
for al in alignments:
    alignments_mp[al["name"]] = al["aln"]

clusters = []
print "Aln num", len(alignments_mp)
for al in alignments_mp:
    new_clusters = []
    first_cl = -1
    in_some_cluster = False
    for c in clusters:
        in_cluster = False
        for c_al in c:
            if aai([alignments_mp[al], alignments_mp[c_al]]) > 90:
                in_cluster = True
                break
        if in_cluster:
            in_some_cluster = True
            if first_cl == -1:
                c.add(al)
                new_clusters.append(c)
                first_cl = len(new_clusters) - 1
            else:
                new_clusters[first_cl] |= c
        else:
            new_clusters.append(c)
    if not in_some_cluster:
        new_clusters.append(set({al}))
    clusters = []
    for cl in new_clusters:
        clusters.append(cl)

print "Number of clusters", len(clusters)
only_amr = 0
only_spaligner = 0
for c in clusters:
    has_amr = False
    has_spaligner = False
    for cc in c:
        has_amr = has_amr or cc.startswith("AMR")
        has_spaligner = has_spaligner or cc.startswith("SP")
    if has_amr and not has_spaligner:
        only_amr += 1
        print "AMR only"
        for cc in c:
            print cc, aai([alignments_mp[cc], proteins[cc[len("AMR_"):]].seq])
            
    if not has_amr and has_spaligner:
        only_spaligner += 1
        print "SP only"
        for cc in c:
            print cc, aai([alignments_mp[cc], proteins[cc[len("SP_"):]].seq])

print "Number of amr clusters", only_amr
print "Number of spaligner clusters", only_spaligner

