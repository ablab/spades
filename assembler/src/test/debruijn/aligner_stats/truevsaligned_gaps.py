import edlib
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

def load_edges(file):
    res = {}
    name = ""
    with open(file, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith(">"):
                name = ln.strip()[1:]
                res[name] = ""
            else:
                res[name] += ln.strip()
    return res

def load_alignments(filename):
    res = {}
    fin = open(filename, "r")
    bwa_num = 0
    for ln in fin.readlines():
        cur_read, seq_start, seq_end, rlen, path_dirty, edgelen = ln.strip().split("\t")
        cur_read = cur_read.split(" ")[0]
        p_lst = []
        for p in path_dirty.split("]")[:-1]:

            lst = p.split(" ")
            if p.startswith(","):
                edgeid = lst[1]
                lst.pop(0)
            else:
                edgeid = lst[0]
            e_start, e_end, r_start, r_end = int(lst[1].split(",")[0][1:]), int(lst[1].split(",")[1][:-1]), int(lst[2].split(",")[0][1:]), int(lst[2].split(",")[1])
            
            p_lst.append({"name": edgeid, "e_start": e_start, "e_end": e_end, "r_start": r_start, "r_end": r_end})
        path = []
        edge_tag = []
        bwa_path = []
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
        res[cur_read] = {"path_with_ends": p_lst, "path": path, "edge_tag": edge_tag, "bwa_path": bwa_path}
    fin.close()
    return res

def load_truealignments(filename, reads):
    res = {}
    fin = open(filename, "r")
    bwa_num = 0
    for ln in fin.readlines():
        cur_read, rlen, pathsz, path_dirty, edgelen, edgelen2, mapped_seq = ln.strip().split("\t")
        cur_read = cur_read.split(" ")[0]
        if cur_read in reads:
            p_lst = []
            for p in path_dirty.split("]")[:-1]:
                lst = p.split(" ")
                if p.startswith(","):
                    edgeid = lst[1]
                    lst.pop(0)
                else:
                    edgeid = lst[0]
                e_start, e_end, r_start, r_end = int(lst[1].split(",")[0][1:]), int(lst[1].split(",")[1][:-1]), int(lst[2].split(",")[0][1:]), int(lst[2].split(",")[1])
                p_lst.append({"name": edgeid, "e_start": e_start, "e_end": e_end, "r_start": r_start, "r_end": r_end})
            path = []
            edge_tag = []
            bwa_path = []
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
            res[cur_read] = {"path_with_ends": p_lst, "path": path, "edge_tag": edge_tag, "bwa_path": bwa_path}
    fin.close()
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

def extract_seq(path, edges):
    res = ""

    for e in path:
        res += edges[e["name"]][e["e_start"]: e["e_end"]]
    return res, path[0]["r_start"], path[-1]["r_end"] 


def count_gaps_ed(trueinfo, info, true_ind, ind, read, edges, is_real_reads):
    truepath = trueinfo["path"]
    path = info["path"]
    strange_ed = 0
    normal = 0
    for j in xrange(2, len(true_ind) - 1):
        if truepath[true_ind[j - 1]: true_ind[j]] != path[ind[j - 1]: ind[j]]:
            if len(path[ind[j - 1]: ind[j]]) > 1:
                true_seq, true_start, true_end = extract_seq(trueinfo["path_with_ends"][true_ind[j - 1]: true_ind[j] + 1], edges) 
                seq, start, end = extract_seq(info["path_with_ends"][ind[j - 1]: ind[j] + 1], edges) 
                #print true_start, true_end
                #print start, end
                if is_real_reads:
                    read_seq = read[start: end]
                else:
                    read_seq = read[true_start: true_end]

                true_ed = edlib.align(read_seq, true_seq, mode="NW" )['editDistance']
                aligned_ed = edlib.align(read_seq, seq, mode="NW" )['editDistance']
                #print "  ED ", true_ed, aligned_ed, len(path[ind[j - 1]: ind[j]])
                if true_ed > aligned_ed:
                    strange_ed += 1
                else:
                    normal += 1
    return strange_ed, normal


def count_gaps_ed_all(truepaths, alignedpaths, edges, reads, is_real_reads):
    unknown = 0
    strange_ed = 0
    normal = 0
    for r in alignedpaths.keys():
        if r not in truepaths:
            continue
        if not is_unique(truepaths[r]["path"], alignedpaths[r]["bwa_path"]):
            unknown += 1
        else:
            true_ind = get_bwa_inds(truepaths[r]["path"], alignedpaths[r]["bwa_path"])
            aligned_ind = get_bwa_inds_aligned(alignedpaths[r]["edge_tag"])
            strange_ed_cur, normal_cur= count_gaps_ed(truepaths[r], alignedpaths[r], true_ind, aligned_ind, reads[r], edges, is_real_reads)
            strange_ed += strange_ed_cur
            normal += normal_cur
    return [strange_ed, normal, unknown]

ideal_reads = load_reads(sys.argv[1])
reads = load_reads(sys.argv[2])
edges = load_edges(sys.argv[3])
truepaths = load_truealignments(sys.argv[4], reads)
alignedpaths = load_alignments(sys.argv[5])

print count_gaps_ed_all(truepaths, alignedpaths, edges, reads, True)
print count_gaps_ed_all(truepaths, alignedpaths, edges, ideal_reads, False)

