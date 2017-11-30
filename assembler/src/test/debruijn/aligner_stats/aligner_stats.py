__author__ = 'tanunia'

import sys
import pysam

def load_reads(filename):
    res = {}
    key = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip().split(" ")[0]
                res[key] = ""
            else:
                res[key] += ln.strip()
    return res


def load_reads_len(filename):
    res = {}
    key = ""
    cur_str = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip().split(" ")[0]
                cur_str = ""
                res[key] = 0
            else:
                cur_str += ln.strip()
                res[key] = len(cur_str)

    return res

def load_good_reads(readsbwamem):
    readssam = pysam.AlignmentFile(readsbwamem, "rb")
    mp = {}
    for r in readssam:
        start = r.reference_start
        end = r.reference_end
        name = r.query_name
        if start == -1 or end is None:
            continue
        if abs((end - start) - reads_len[name]) * 1.0 / reads_len[name] < 0.2:
            if name not in mp or abs((end - start) - reads_len[name]) < abs(mp[name]["len"] - reads_len[name]):
                if name not in mp:
                    mp[name] = {}
                mp[name]["len"] = end - start
                mp[name]["start"] = start
                mp[name]["end"] = end
    print "Map loaded #items: ", len(mp)
    names_set = set(mp.keys())
    return names_set

def is_aligned(al_len, str_len):
    if abs(al_len - str_len) * 1.0 / str_len >= 0.4:
        return False
    else:
        return True


if (len(sys.argv) < 4):
    print "Usage: aligner_stats.py <file with .bam> <file with reads .fasta> <file with alignment info>"
    exit(-1)

align_file = sys.argv[3]
read_file = sys.argv[2]
readsbwamem = sys.argv[1]

reads_len = load_reads_len(read_file)
names_set = load_good_reads(readsbwamem)
print "Reads loaded"

total_readlen = 0
total_readnum = 0
for k in reads_len.keys():
    if k in names_set:
        total_readlen += reads_len[k]
        total_readnum += 1

print "Count total len ", total_readlen, total_readnum

aligned_len = 0
aligned_num = 0
path_cnt = 0
path_cnt_len = 0
path_med_len = []
edge_med_len = []

with open(align_file, "r") as fin:
    ln = fin.readline()
    while ln != "":
        [name, start, end, sz, path, path_len, score, subread] = ln.split("\t")
        start = int(start)
        end = int(end)
        name = name.split(" ")[0]
        if name in names_set and is_aligned(end-start, reads_len[name]):
            aligned_num += 1
            aligned_len += (end - start)
            if len(path_len.split(",")) > 2:
                path_cnt += 1
                path_cnt_len += sum([int(x) for x in path_len.split(",")[:-1]])
                path_med_len.append(len(path_len.split(",")) - 1)
                edge_med_len.extend([int(x) for x in path_len.split(",")[:-1]])
        ln = fin.readline()

print "aligned: ", aligned_num, " total: ", total_readnum, " %: ", aligned_num*1.0/total_readnum*100
print "aligned len: ", aligned_len, " total len: ", total_readlen, " %: ", aligned_len*1.0/total_readlen*100
print "paths with more than 1 edge: ", path_cnt, " total len: ", path_cnt_len, " %: ", path_cnt*1.0/total_readnum*100
print "median path len: ", sorted(path_med_len)[len(path_med_len)/2]
print "median edge len: ", sorted(edge_med_len)[len(edge_med_len)/2]
print "avg path len: ", sum(path_med_len)/len(path_med_len)
print "avg edge len: ", sum(edge_med_len)/len(edge_med_len)
