import sys
import re

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

logname = sys.argv[1]
idealfile = sys.argv[2]
thresholds = [0, 10, 20, 50, 100, 200, 300, 400, 500, 700, 900, 1000, 1500]

def load_reads(filename):
    res = {}
    cur_read = ""
    fin = open(filename, "r")
    for ln in fin.readlines():
        if ln.startswith(">"):
            cur_read = ln.strip()[1:].split(" ")[0]
            res[cur_read] = ""
        else:
            res[cur_read] += ln.strip()
    fin.close()
    return res

def load_la_paths(logname):
    edge_sets = {}
    with open(logname, "r") as fin:
        s = fin.readline()
        pattern = re.compile("ReadName=S\d_\d+ BWA set edge=\d+ edge_length=\d+ e_start=\d+ e_end=\d+ r_start=\d+ r_end=\d+")
        while s != "":
            res = re.findall(pattern, s.strip())
            for it in res:
                lst = it.split(" ")
                readname = lst[0][len("ReadName="):]
                edgeid = lst[3][len("edge="):]
                edge_len = int(lst[4][len("edge_length="):])
                e_start = int(lst[5][len("e_start="):])
                e_end = int(lst[6][len("e_end="):])
                r_start = int(lst[7][len("r_start="):])
                r_end = int(lst[8][len("r_end="):])
                if readname not in edge_sets:
                    edge_sets[readname] = set()
                edge_sets[readname].add(edgeid) #.append([edgeid, str(e_end - e_start), 2*(e_end - e_start) < min(r_start, e_start) + min(edge_len - e_end, reads_len[readname] - r_end)])
                #print [edgeid, str(e_end - e_start), ]
            s = fin.readline()

    return edge_sets

def load_ideal_paths(idealfile):
    read_set, read_set_len = {}, {}
    with open(idealfile, "r") as fin:
        for ln in fin.readlines():
            [name, nuc_len, path_len, path, edge_len, whole_edge_len, aln] = ln.strip().split("\t")
            path_lst = []
            for x in path.split("]")[:-1]:
                if x.startswith(","):
                    path_lst.append(x.split()[1])
                else:
                    path_lst.append(x.split()[0])
            edge_len_lst = edge_len.split(",")[:-1]
            whole_edge_len_lst = [int(x) for x in whole_edge_len.split(",")[:-1]]
            read_set[name] = path_lst
            read_set_len[name] = whole_edge_len_lst

    return read_set, read_set_len

def good_length(l, ths):
    if l >= ths[0] and l < ths[1]:
        return True
    else:
        return False

def filter_ideal_paths(paths, lens, ths):
    edge_sets = {}
    for name in paths.keys():
        whole_edge_len_lst = lens[name]
        path_lst = paths[name]
        add = False
        for e in whole_edge_len_lst:
            if good_length(e, ths):
                add = True
                break
        if name not in edge_sets:
            edge_sets[name] = set()
        if add:
            for i in xrange(len(path_lst)):
                if good_length(whole_edge_len_lst[i], ths):
                    edge_sets[name].add(path_lst[i])
    return edge_sets

reads = load_reads(sys.argv[3])
la_paths = load_la_paths(logname)
ideal_paths_all, ideal_paths_lens = load_ideal_paths(idealfile)
for ind in xrange(1, len(thresholds)):
    threshold1, threshold2 = thresholds[ind - 1], thresholds[ind]
    ideal_paths = filter_ideal_paths(ideal_paths_all, ideal_paths_lens, [threshold1, threshold2])
    fn = 0
    total = 0
    for k in reads.keys():
        if k not in ideal_paths:
            continue
        total += len(ideal_paths[k])
        if k not in la_paths:
            fn += len(ideal_paths[k])
        else:
            #print ideal_paths[k]
            #print la_paths[k]
            fn += len(ideal_paths[k] - la_paths[k])

    print "Length between ", threshold1, "-", threshold2, " missing hits", (fn*1.0/total)*100, "% (", fn, "out of", total, ")"
        