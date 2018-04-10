import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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

def load_truepaths(filename):
    res = {}
    fin = open(filename, "r")
    for ln in fin.readlines():
        cur_read, rlen, pathsz, path_dirty, edgelen, edgelen2, mapped_seq = ln.strip().split("\t")
        cur_read = cur_read.split(" ")[0]
        path = []
        seq = []
        for x in path_dirty.split("]")[:-1]:
            if x.startswith(","):
                path.append(x.split()[1])
            else:
                path.append(x.split()[0])
            seq.append(x.split("[")[1].split(",")[0])

        res[cur_read] = { "len": int(rlen), "path": path, "edgelen": [int(x) for x in edgelen.split(",")[:-1]], "mapped": mapped_seq, "seq": seq }
    fin.close()
    return res


def load_alignments(filename):
    res = {}
    fin = open(filename, "r")
    bwa_num = []
    for ln in fin.readlines():
        #cur_read, seq_start, seq_end, rlen, path_dirty, edgelen, ed = ln.strip().split("\t")
        cur_read, seq_start, seq_end, rlen, path_dirty, edgelen, ss, ed = ln.strip().split("\t")
        cur_read = cur_read.split(" ")[0]
        path = []
        edge_tag = []
        bwa_path = []
        seq = []
        seq_ends = []
        bwa_num.append(0)
        for x in path_dirty.split("]")[:-1]:
            if x.startswith(","):
                path.append(x.split()[1])
            else:
                path.append(x.split()[0])
            if "[0,0" in x:
                edge_tag.append("in")
            else:
                edge_tag.append("bwa")
                bwa_num[len(bwa_num) - 1] += 1
                if x.startswith(","):
                    bwa_path.append(x.split()[1])
                else:
                    bwa_path.append(x.split()[0])
            seq.append(x.split("[")[1].split(",")[0])
            seq_ends.append(x.split("[")[1].split(",")[1])
        #if int(seq_end) - int(seq_start) > 1000:
        res[cur_read] = { "len": int(rlen), "path": path, "bwa_path": bwa_path, "edgelen": edgelen.split(",")[:-1], "edge_tag": edge_tag, \
                            "mapped_s":int(seq_start), "mapped_e":int(seq_end), "seq": seq , "seq_end": seq_ends }
    fin.close()
    print "Median number of bwa hits per read:", sorted(bwa_num)[len(bwa_num)/2]
    return res


def cnt_badideal(reads, alignedpaths):
    res = 0
    for r in reads.keys():
        if r not in alignedpaths:
            res += 1
    return res

def cnt_notmapped(reads, truepaths, alignedpaths):
    res = 0
    for r in reads.keys():
        if r not in truepaths:
            continue
        if r not in alignedpaths:
            #print r
            res += 1
    return res

def cnt_problempaths(truepaths, alignedpaths):
    res = 0
    for r in alignedpaths.keys():
        if r not in truepaths.keys():
            continue
        if truepaths[r]["path"] != alignedpaths[r]["path"]:
            res += 1
    return res

def cnt_problembwa(truepaths, alignedpaths):
    res = 0
    total_path = 0
    for r in alignedpaths.keys():
        if r not in truepaths:
            continue
        j = 0
        total_path += 1
        for e in truepaths[r]["path"]:
            if j < len(alignedpaths[r]["bwa_path"]) and e == alignedpaths[r]["bwa_path"][j]:
                j += 1
        if j < len(alignedpaths[r]["bwa_path"]):
            res += 1
    return res

def get_bwa_inds(path, bwapath):
    ind = []
    j = 0
    for i in xrange(len(path)):
        if j < len(bwapath) and path[i] == bwapath[j]:
            j += 1
            ind.append(i)
    ind.insert(0, 0)
    ind.append(len(path))
    return ind

def get_bwa_inds_aligned(path):
    ind = []
    for i in xrange(len(path)):
        if path[i] == "bwa":
            ind.append(i)
    ind.insert(0, 0)
    ind.append(len(path))
    return ind


def is_unique(path, subpath):
    ind = get_bwa_inds(path, subpath)
    if len(ind) < len(subpath) + 2:
        return False
    for j in xrange(1, len(subpath) + 1):
        if subpath[j - 1] in path[ind[j - 1]: ind[j]] or subpath[j - 1] in path[ind[j] + 1: ind[j + 1]]:
            return False
    return True

def is_wrong_start(truepath, path, true_ind, ind, edgelen):
    empty = False
    if truepath[true_ind[0]: true_ind[1]] == path[ind[0]: ind[1]] :
        return [False, empty, -1]
    else:
        if len(path[ind[0]: ind[1]]) == 0:
            empty = True
        ln = 0
        for it in xrange(true_ind[0], true_ind[1]):
            ln += edgelen[it]
        #print truepath[true_ind[0]: true_ind[1]], path[ind[0]: ind[1]] 
        return [True, empty, ln]

def is_wrong_end(truepath, path, true_ind, ind, edgelen):
    empty = False
    if truepath[true_ind[-2]: true_ind[-1]] == path[ind[-2]: ind[-1]] :
        return [False, empty, -1]
    else:
        if len(path[ind[-2]: ind[-1]]) == 1:
            empty = True
        ln = 0
        for it in xrange(true_ind[-2], true_ind[-1]):
            ln += edgelen[it]
        return [True, empty, ln]


def form_length_lists(tpaths, apaths, K, edges = True, use_edges_len =True):
    unknown = 0
    res_prefix = []
    res_suffix = []
    res_plus = []
    wrong_start = 0
    wrong_end = 0
    for r in apaths.keys():
        if r not in tpaths:
            continue
        if not is_unique(tpaths[r]["path"], apaths[r]["bwa_path"]):
            unknown += 1
        else:
            if edges:
                true_ind = get_bwa_inds(tpaths[r]["path"], apaths[r]["bwa_path"])
                aligned_ind = get_bwa_inds_aligned(apaths[r]["edge_tag"])
                has_wrong_start, empty, ln = is_wrong_start(tpaths[r]["path"], apaths[r]["path"], true_ind, aligned_ind, tpaths[r]["edgelen"])
                m = 0
                if has_wrong_start:
                    wrong_start += 1
                    if use_edges_len:
                        res_prefix.append(ln)
                        m += ln
                    else:    
                        res_prefix.append( apaths[r]["mapped_s"] )
                        m += apaths[r]["mapped_s"] 
                has_wrong_end, empty, ln = is_wrong_end(tpaths[r]["path"], apaths[r]["path"], true_ind, aligned_ind, tpaths[r]["edgelen"])
                if has_wrong_end:
                    wrong_end += 1
                    if use_edges_len:
                        res_suffix.append(ln)
                        m += ln
                    else:
                        res_suffix.append( apaths[r]["len"] - apaths[r]["mapped_e"] - K)
                        m += apaths[r]["len"] - apaths[r]["mapped_e"] - K
                if has_wrong_end or has_wrong_start:
                    res_plus.append(m)
            else:
                if apaths[r]["mapped_s"] > 0:
                    res_prefix.append( apaths[r]["mapped_s"] )
                if apaths[r]["len"] - apaths[r]["mapped_e"] - K > 0:
                    res_suffix.append( apaths[r]["len"] - apaths[r]["mapped_e"] - K)
            #if apaths[r]["len"] - apaths[r]["mapped_e"] - K + apaths[r]["mapped_s"] > 100:
                if apaths[r]["mapped_s"] > 0 or apaths[r]["len"] - apaths[r]["mapped_e"] - K > 0:
                    res_plus.append( apaths[r]["len"] - apaths[r]["mapped_e"] - K + apaths[r]["mapped_s"])

    print "Unknown=", unknown
    print sum(res_prefix)*1.0/len(res_prefix), sum(res_suffix)*1.0/len(res_suffix), sum(res_plus)*1.0/len(res_plus)
    print sorted(res_prefix)[len(res_prefix)/2], sorted(res_suffix)[len(res_suffix)/2], sorted(res_plus)[len(res_plus)/2] 
    return res_prefix, res_suffix, res_plus



def build_histogram(tpaths, apaths1, apaths2, K, name):
    res_prefix1, res_suffix1, res_plus1 = form_length_lists(tpaths, apaths1, K)
    res_prefix2, res_suffix2, res_plus2 = form_length_lists(tpaths, apaths2, K)
    plt.figure()
    plt.hist(res_plus1, 80, color="blue", alpha=0.4, label="initial")
    plt.hist(res_plus2, 80, color="green", alpha=0.4, label="new weights")
    #plt.axis([0, 1000, 0, 15000])
    plt.title("Lost on edges: " + name)
    plt.legend(loc='upper right')
    plt.savefig(name + '_edge_length_hist_sum.png')

    plt.figure()
    plt.hist(res_prefix1, 80, color="blue", alpha=0.4, label="initial")
    plt.hist(res_prefix2, 80, color="green", alpha=0.4, label="new weights")
    #plt.axis([0, 1000, 0, 1000])
    plt.title("Lost in prefix: " + name)
    plt.legend(loc='upper right')
    plt.savefig(name + '_edge_length_hist_prefix.png')

    plt.figure()
    plt.hist(res_suffix1, 80, color="blue", alpha=0.4, label="initial")
    plt.hist(res_suffix2, 80, color="green", alpha=0.4, label="new weights")
    #plt.axis([0, 1000, 0, 1000])
    plt.title("Lost in suffix: " + name)
    plt.legend(loc='upper right')
    plt.savefig(name + '_edge_length_hist_suffix.png')

reads = load_reads(sys.argv[1])
truepaths = load_truepaths(sys.argv[2])
alignedpaths1 = load_alignments(sys.argv[3])
alignedpaths2 = load_alignments(sys.argv[4])
K = int(sys.argv[5])
name = sys.argv[6]

print "Total=", len(reads), " ideal=", len(reads) - cnt_badideal(reads, truepaths)
print "Notmapped1=", cnt_notmapped(reads, truepaths, alignedpaths1), " Notmapped2=", cnt_notmapped(reads, truepaths, alignedpaths2) 
print "Warning: works only without GrowEnds()"
print "Draw histrogram of distribution of lost edges length.."
build_histogram(truepaths, alignedpaths1, alignedpaths2, K, name) 








