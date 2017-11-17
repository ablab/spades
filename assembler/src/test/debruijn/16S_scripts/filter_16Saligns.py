
__author__ = 'tanunia'

import sys
import re
import multiprocessing

import edlib
import editdistance


def load_reads(filename):
    res = {}
    key = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip()
                res[key] = ""
            else:
                res[key] += ln.strip()
    return res


def is_aligned(al_len, str_len):
    if abs(al_len - str_len) * 1.0 / str_len >= 0.4:
        return False
    else:
        return True

def edist(lst):
    #return editdistance.eval(lst[0], lst[1])
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

def load_aligns(alnfile, reads):
    paths = []
    with open(alnfile, "r") as aln:
        for ln in aln.readlines():
            cur_read, seq_start, seq_end, read_len, path, edges_len, score, align = ln.strip().split("\t")
            if is_aligned(len(align), int(read_len)):
                paths.append({"path": path, "start": int(seq_start), "end": int(seq_end),
                              "edge_lens": edges_len.split(",")[:-1], "read": cur_read, "read_len": int(read_len),
                              "score": int(score)*1.0/int(read_len), "score2": int(score), "aln": align})

    # pool = multiprocessing.Pool(16)
    # #scores = pool.map(edist, zip(["bahama"], ["banana"]))
    # scores = pool.map(edist, zip([reads[path["read"]] for path in paths], [path["aln"] for path in paths]))
    # for i in xrange(len(paths)):
    #     paths[i]["score"] = scores[i]
    return paths


def load_names():
    res = set()
    with open("/Sid/tdvorkina/gralign/16S/synth16S_names.txt", "r") as inf:
        for name in inf.readlines():
            res.add(name.strip())
    print "Loaded ", len(res), " orgs"
    return res


if len(sys.argv) < 3:
    print "Usage: filter_16Saligns.py <file with reads .fasta> <file with alignment info> <out-fasta>"
    exit(-1)

readsfile= sys.argv[1]
alnfile = sys.argv[2]
outfile = sys.argv[3]
names = load_names()

reads = load_reads(readsfile)
print "Reads loaded"

paths = load_aligns(alnfile, reads)
print "Paths loaded"

paths = sorted(paths,  key=lambda x: x["score"])
cnt = 0
used = set()
lens = []
plens = []
scores = []
orgs = set()
with open(outfile, "w") as fout:
    for p in paths:
        to_add = False
        cur_name = p["read"]
        scores.append(p["score"])
        if p["score"] < 0.2:
            cur_path = re.findall("^\d+| \d+", p["path"])
            sm = 0
            unique_edge = ""
            unique_edge_len = ""
            for ind in xrange(len(cur_path)):
                e = cur_path[ind]
                e_len = int(p["edge_lens"][ind])
                if e not in used and int(e_len) >= 300:
                    to_add = True
                    unique_edge = e
                    unique_edge_len = int(e_len)
                    break
            to_add = True
            if to_add:
                nm = ''
                # for name in names:
                #     if p["read"].startswith(name):
                #         orgs.add(name)
                #         nm = name
                for ind in xrange(len(cur_path)):
                    e = cur_path[ind]
                    used.add(e)
                    e_len = int(p["edge_lens"][ind])
                    lens.append(e_len)
                    plens.append(len(cur_path))
                cnt += 1
                #print p["read"], p["score2"], p["score"], len(p["aln"]), p["read_len"], ",".join(cur_path)
                fout.write(">" + p["read"] + "\n" + p["aln"] + "\n")
            # else:
            #     print "FAILED ", p["read"], p["score"], len(p["aln"]), p["read_len"]

print "\n".join(sorted(orgs))
print "16S=", cnt, " Paths=", len(paths), " Reads=", len(reads)
print "Orgs num=", len(orgs),
print "Median edge len: ", sorted(lens)[len(lens)/2]
print "Median path len(in edges): ", sorted(plens)[len(plens)/2]
print "Median alignment score: ",sorted(scores)[len(scores)/2]
#print "\n".join(sorted(names - orgs))
