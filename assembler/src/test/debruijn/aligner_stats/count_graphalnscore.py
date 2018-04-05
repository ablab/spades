__author__ = 'tanunia'

import sys
import pysam
import subprocess
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO
import pandas as pd


def load_reads(filename):
    res = {}
    key = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip().split(" ")[0]
                res[key] = 0
            else:
                res[key] += len(ln.strip())
    return res

def form_fasta_forbwamem(alnfile, prefix):
    alns = {}
    reads_set = set()
    with open(prefix + ".fasta", "w+") as forbwamem:
        with open(alnfile, "r") as alnf:
            for ln in alnf.readlines():
                lst = ln.split("\t")
                read_name = lst[0].split(" ")[0]
                aln = lst[-2].strip()
                alns[read_name] = len(aln)
                reads_set.add(read_name)
                forbwamem.write(">" + read_name + "\n" + aln + "\n")
    return prefix + ".fasta", alns

def load_good_reads(readsbwamem):
    readssam = pysam.AlignmentFile(readsbwamem, "rb")
    mp = {}
    for r in readssam:
        start = r.reference_start
        end = r.reference_end
        name = r.query_name
        if start == -1 or end is None:
            continue
        if abs((end - start) - reads[name]) * 1.0 / reads[name] < 0.2:
            if name not in mp or abs((end - start) - reads[name]) < abs(mp[name]["len"] - reads[name]):
                if name not in mp:
                    mp[name] = {}
                mp[name]["len"] = end - start
                mp[name]["start"] = start
                mp[name]["end"] = end
    print "Map loaded #items: ", len(mp)
    return mp

def is_aligned(al_len, str_len):
    if abs(al_len - str_len) * 1.0 / str_len >= 0.4:
        return False
    else:
        return True




if len(sys.argv) < 5:
    print "Usage: count_graphalnscore.py <output-prefix> <graph_aln_file> <ref> <reads> {reads-bam}"
    exit(-1)

print "Form fasta from alignments"
prefix = sys.argv[1]
alnfile = sys.argv[2]
forbwamemf, alns = form_fasta_forbwamem(alnfile, prefix)

print "Load reads"
readsfile = sys.argv[4]
reads = load_reads(readsfile)

print "Run bwa mem on alignments"
ref = sys.argv[3]
alnsbwamem = prefix + "_aln.bam"
subprocess.call(["bwa mem -x pacbio -t 33 " + ref + " " + forbwamemf + " | samtools view -bS - > " + alnsbwamem], shell=True)

readsbwamem = ""
if len(sys.argv) == 5:
    print "Run bwa mem on origin reads"
    readsbwamem = prefix + "_origin.bam"
    subprocess.call(["bwa mem -x pacbio -t 33 " + ref + " " + readsfile + " | samtools view -bS - > " + readsbwamem], shell=True)
else:
    readsbwamem = sys.argv[5]

print "Load perfectly mapped reads"
mp = load_good_reads(readsbwamem)

print "Filter good reads alignments"
alnssam = pysam.AlignmentFile(alnsbwamem, "rb")

buckets = {0.05: 0, 0.1: 0, 0.2: 0, 0.3: 0, 0.4: 0, 0.5: 0, 0.6: 0, 0.7:0, 0.8: 0, 0.9: 0, 1: 0, 1.5: 0, 100500000: 0}
total = set()
total_aligned_with_ga = set()
mp_aligned = {}
for aln in alnssam:
    start = aln.reference_start
    end = aln.reference_end
    name = aln.query_name
    if start == -1 or end == None:
        continue
    if name in mp:
        total.add(name)
        if is_aligned(alns[name], reads[name]):
            total_aligned_with_ga.add(name)
            if name not in mp_aligned or abs((end - start) - reads[name]) < abs(mp_aligned[name]["len"] - reads[name]):
                if name not in mp_aligned:
                    mp_aligned[name] = {}
                for k in sorted(buckets.keys()):
                    if abs(end - mp[name]["end"]) < k*reads[name] and abs(start - mp[name]["start"]) < k*reads[name]:
                        mp_aligned[name]["len"] = (end - start)
                        mp_aligned[name]["bucket"] = k
                        break

for k in mp_aligned.keys():
    buckets[mp_aligned[k]["bucket"]] += 1

print "Total reads number: ", len(reads), "; Aligned on ref reads number: ", len(mp)
print "Total graphaln number: ", len(total)*1.0/len(mp), "; Total aligned with ga: ",len(total_aligned_with_ga)*1.0/len(mp) , "; Filtered by aligned reads: ", len(mp_aligned)*1.0/len(mp)
sm = 0
for k in sorted(buckets.keys()):
    sm += buckets[k]
    print k, sm, sm*1.0/len(mp)
