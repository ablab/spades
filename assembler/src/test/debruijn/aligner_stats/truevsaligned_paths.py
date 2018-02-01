import sys


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

        res[cur_read] = { "len": rlen, "path": path, "edgelen": edgelen.split(",")[:-1], "mapped": mapped_seq, "seq": seq }
    fin.close()
    return res


def load_alignments(filename):
    res = {}
    fin = open(filename, "r")
    bwa_num = 0
    for ln in fin.readlines():
        cur_read, seq_start, seq_end, rlen, path_dirty, edgelen, trash1, trash2 = ln.strip().split("\t")
        cur_read = cur_read.split(" ")[0]
        path = []
        edge_tag = []
        bwa_path = []
        seq = []
        seq_ends = []
        for x in path_dirty.split("]")[:-1]:
            if x.startswith(","):
                path.append(x.split()[1])
            else:
                path.append(x.split()[0])
            if "[0,0" in x:
                edge_tag.append("in")
            else:
                edge_tag.append("bwa")
                bwa_num += 1
                if x.startswith(","):
                    bwa_path.append(x.split()[1])
                else:
                    bwa_path.append(x.split()[0])
            seq.append(x.split("[")[1].split(",")[0])
            seq_ends.append(x.split("[")[1].split(",")[1])

        res[cur_read] = { "len": rlen, "path": path, "bwa_path": bwa_path, "edgelen": edgelen.split(",")[:-1], "edge_tag": edge_tag, \
                            "mapped_s":int(seq_start), "mapped_e":int(seq_end), "seq": seq , "seq_end": seq_ends }
    fin.close()
    print "Number of edges detected by bwa:", bwa_num
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

def cover_same_subseq(a_start, a_end, b_start, b_end):
    OVERLAP = 20
    if ( ((b_end - a_start)*(b_start - a_start) < 0 and (b_end - a_start) < OVERLAP) or (b_end - a_start)*(b_start - a_start) > 0) \
             and  ( ((b_end - a_end)*(b_start - a_end) < 0 and (a_end - b_start) < OVERLAP) or (b_end - a_end)*(b_start - a_end) > 0) \
             and  ( ((a_end - b_start)*(a_start - b_start) < 0 and (a_end - b_start) < OVERLAP) or (a_end - b_start)*(a_start - b_start) > 0) \
             and  ( ((a_end - b_end)*(a_start - b_end) < 0 and (b_end - a_start) < OVERLAP) or (a_end - b_end)*(a_start - b_end) > 0)  :
        return false;
    return true;

def cnt_problembwa(truepaths, alignedpaths):
    res = 0
    for r in alignedpaths.keys():
        if r not in truepaths:
            continue
        j = 0
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


def is_wrong_start(truepath, path, true_ind, ind):
    empty = False
    if truepath[true_ind[0]: true_ind[1]] == path[ind[0]: ind[1]] :
        return [False, empty]
    else:
        if len(path[ind[0]: ind[1]]) == 0:
            empty = True
        return [True, empty]

def is_wrong_end(truepath, path, true_ind, ind):
    empty = False
    if truepath[true_ind[-2]: true_ind[-1]] == path[ind[-2]: ind[-1]] :
        return [False, empty]
    else:
        if len(path[ind[-2]: ind[-1]]) == 1:
            empty = True
        return [True, empty]

def is_wrong_gap(truepath, path, true_ind, ind):
    empty = False
    has_wrong_gap = False
    for j in xrange(2, len(true_ind) - 1):
        if truepath[true_ind[j - 1]: true_ind[j]] != path[ind[j - 1]: ind[j]]:
            if len(path[ind[j - 1]: ind[j]]) == 1:
                empty = True
            else:
                has_wrong_gap = True
    return [has_wrong_gap, empty]


def cnt_wronglyclosedgaps(truepaths, alignedpaths):
    res, res_start, res_end = 0, 0, 0
    res_empty, res_start_empty, res_end_empty = 0, 0, 0
    unknown = 0
    for r in alignedpaths.keys():
        if r not in truepaths:
            continue
        if not is_unique(truepaths[r]["path"], alignedpaths[r]["bwa_path"]):
            unknown += 1
        else:
            true_ind = get_bwa_inds(truepaths[r]["path"], alignedpaths[r]["bwa_path"])
            aligned_ind = get_bwa_inds_aligned(alignedpaths[r]["edge_tag"])
            has_wrong_start, empty = is_wrong_start(truepaths[r]["path"], alignedpaths[r]["path"], true_ind, aligned_ind)
            if has_wrong_start:
                if empty:
                    res_start_empty += 1
                else:
                    res_start += 1

            has_wrong_end, empty = is_wrong_end(truepaths[r]["path"], alignedpaths[r]["path"], true_ind, aligned_ind)
            if has_wrong_end:
                if empty:
                    res_end_empty += 1
                else:
                    res_end += 1
            has_wrong_gap, empty = is_wrong_gap(truepaths[r]["path"], alignedpaths[r]["path"], true_ind, aligned_ind)
            if has_wrong_gap:
                res += 1
            if empty:
                res_empty += 1
    return [res, res_empty, res_start, res_start_empty, res_end, res_end_empty, unknown]

reads = load_reads(sys.argv[1])
truepaths = load_truepaths(sys.argv[2])
alignedpaths = load_alignments(sys.argv[3])

print "Total=", len(reads), " notideal=", cnt_badideal(reads, truepaths),  " notmapped=", cnt_notmapped(reads, truepaths, alignedpaths) 
print "Paths with problems ", cnt_problempaths(truepaths, alignedpaths)
print "BWA fail ", cnt_problembwa(truepaths, alignedpaths)
print "Gaps stage problems (in, prefix, suffix, unknown) ", cnt_wronglyclosedgaps(truepaths, alignedpaths) 








