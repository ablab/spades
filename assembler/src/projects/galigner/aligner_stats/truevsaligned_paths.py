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

        res[cur_read] = { "len": rlen, "path": path, "edgelen": [int(x) for x in edgelen.split(",")[:-1]], "mapped": mapped_seq, "seq": seq }
    fin.close()
    return res


def load_alignments(filename):
    res = {}
    fin = open(filename, "r")
    bwa_num = 0
    for ln in fin.readlines():
        cur_read, seq_start, seq_end, rlen, path_dirty, edgelen, ed = ln.strip().split("\t")
        #cur_read, seq_start, seq_end, rlen, path_dirty, edgelen, ss, ed = ln.strip().split("\t")
        cur_read = cur_read.split(" ")[0]
        path = []
        edge_tag = []
        bwa_path = []
        seq = []
        seq_ends = []
        empty = path_dirty.count(";")
        path_dirty = path_dirty.replace(";", "")
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
        #if int(seq_end) - int(seq_start) > 1000:
        res[cur_read] = { "len": rlen, "path": path, "bwa_path": bwa_path, "edgelen": edgelen.split(",")[:-1], "edge_tag": edge_tag, \
                            "mapped_s":int(seq_start), "mapped_e":int(seq_end), "seq": seq , "seq_end": seq_ends, "empty": empty }
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
            # print r
            # print truepaths[r]["path"]
            # print alignedpaths[r]["bwa_path"] 
            # print alignedpaths[r]["path"]
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


def is_wrong_start(truepath, path, true_ind, ind, edgelen, K):
    empty = False
    if truepath[true_ind[0]: true_ind[1]] == path[ind[0]: ind[1]] :
        return [False, empty]
    else:
        # s = 0
        # for j in xrange(true_ind[0], true_ind[1] + 1):
        #     if truepath[j: true_ind[1]] == path[ind[0]: ind[1]] and s < K + 1:
        #         return [False, empty]
        #     if j < len(edgelen):
        #         s += edgelen[j]
        if len(path[ind[0]: ind[1]]) == 0:
            empty = True
        return [True, empty]

def is_wrong_end(truepath, path, true_ind, ind, edgelen, K):
    empty = False
    if truepath[true_ind[-2]: true_ind[-1]] == path[ind[-2]: ind[-1]] :
        return [False, empty]
    else:
        # for j in xrange(true_ind[-2] + 1, true_ind[-1] + 1):
        #     s = 0
        #     for it in edgelen[j:]:
        #         s += it
        #     if truepath[true_ind[-2]: j] == path[ind[-2]: ind[-1]] and s < K + 1:
        #         print "Found wrong end"
        #         print truepath[true_ind[-2]: true_ind[-1]]
        #         print path[ind[-2]: ind[-1]] 
        #         return [False, empty]

        if len(path[ind[-2]: ind[-1]]) == 1:
            empty = True
        return [True, empty]

def is_wrong_gap(truepath, path, empty, true_ind, ind):
    wrong_gap = 0
    for j in xrange(2, len(true_ind) - 1):
        if truepath[true_ind[j - 1]: true_ind[j]] != path[ind[j - 1]: ind[j]]:
            wrong_gap += 1
    return wrong_gap - empty > 0


def is_wrong_path(truepath, path):
    has_good_path = False
    for j in xrange(len(truepath) - len(path) + 1):
        if " ".join(truepath[j: j + len(path)]) == " ".join(path):
            has_good_path = True
    return has_good_path

def cnt_wronglyclosedgaps(truepaths, alignedpaths, K):
    res, res_start, res_end = 0, 0, 0
    res_empty, res_start_empty, res_end_empty, res_good = 0, 0, 0, 0
    unknown = 0
    for r in alignedpaths.keys():
        if r not in truepaths:
            continue
        if not is_unique(truepaths[r]["path"], alignedpaths[r]["bwa_path"]):
            unknown += 1
        else:
            true_ind = get_bwa_inds(truepaths[r]["path"], alignedpaths[r]["bwa_path"])
            aligned_ind = get_bwa_inds_aligned(alignedpaths[r]["edge_tag"])
            has_wrong_start, empty = is_wrong_start(truepaths[r]["path"], alignedpaths[r]["path"], true_ind, aligned_ind, truepaths[r]["edgelen"], K)
            if has_wrong_start:
                if empty:
                    res_start_empty += 1
                else:
                    res_start += 1

            has_wrong_end, empty = is_wrong_end(truepaths[r]["path"], alignedpaths[r]["path"], true_ind, aligned_ind, truepaths[r]["edgelen"], K)
            if has_wrong_end:
                if empty:
                    res_end_empty += 1
                else:
                    res_end += 1
            has_wrong_gap = is_wrong_gap(truepaths[r]["path"], alignedpaths[r]["path"], alignedpaths[r]["empty"], true_ind, aligned_ind)
            if has_wrong_gap:
                res += 1
            if alignedpaths[r]["empty"] > 0:
                res_empty += 1
    return [res, res_empty, res_start, res_start_empty, res_end, res_end_empty, unknown]

def is_wrong_start2(truepath, path, true_ind, ind, edgelen):
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

def is_wrong_end2(truepath, path, true_ind, ind, edgelen):
    empty = False
    if truepath[true_ind[-2]: true_ind[-1]] == path[ind[-2]: ind[-1]] :
        return [False, empty, -1]
    else:
        if len(path[ind[-2]: ind[-1]]) == 1:
            empty = True
        ln = 0
        for it in xrange(true_ind[-2] + 1, true_ind[-1]):
            ln += edgelen[it]
        return [True, empty, ln]


def cnt_median_alignment_length(tpaths, apaths, K, edges = True, use_edges_len =True):
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
                has_wrong_start, empty, ln = is_wrong_start2(tpaths[r]["path"], apaths[r]["path"], true_ind, aligned_ind, tpaths[r]["edgelen"])
                m = 0
                if has_wrong_start:
                    wrong_start += 1
                    if use_edges_len:
                        res_prefix.append(ln)
                        m += ln
                    else:    
                        res_prefix.append( apaths[r]["mapped_s"] )
                        m += apaths[r]["mapped_s"] 
                has_wrong_end, empty, ln = is_wrong_end2(tpaths[r]["path"], apaths[r]["path"], true_ind, aligned_ind, tpaths[r]["edgelen"])
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

    # print "Unknown=", unknown
    # print sum(res_prefix)*1.0/len(res_prefix), sum(res_suffix)*1.0/len(res_suffix), sum(res_plus)*1.0/len(res_plus)
    return sorted(res_prefix)[len(res_prefix)/2], sorted(res_suffix)[len(res_suffix)/2], sorted(res_plus)[len(res_plus)/2] 

def make_table(results, row_names, caption, name):
    html = """<html><table border="1"><caption>{}</caption><tr><th></th>""".format(name)
    for run_name in sorted(results.keys()):
        html += """<th><div style="width: 200px; height: 50px; overflow: auto">{}</div></th>""".format(run_name)
    html += "</tr>"
    for stat in row_names:
        html += "<tr><td>{}</td>".format(stat)
        for run_name in sorted(results.keys()):
            item = results[run_name]
            html += "<td>{}</td>".format(str(item[stat]))
        html += "</tr>"
    html += "</table>"
    html += "<p>{}</p>".format("<br>".join(caption))
    html += "</html>"
    return html

def save_html(s, fl):
    with open(fl, "w") as fout:
        fout.write(s)

reads = load_reads(sys.argv[1])
truepaths = load_truepaths(sys.argv[2])
aligned_files = sys.argv[3].strip().split(",")
K = int(sys.argv[4])
html_name = sys.argv[5]
aligned_onref = load_reads(sys.argv[6])
res = {}
for fl in aligned_files:
    alignedpaths = load_alignments(fl)
    print "Total=", len(reads), " ideal=", len(reads) - cnt_badideal(reads, truepaths),  " notmapped=", cnt_notmapped(reads, truepaths, alignedpaths) 
    path_problems = cnt_problempaths(truepaths, alignedpaths)
    print "Paths with problems ", path_problems
    bwa_problems = cnt_problembwa(truepaths, alignedpaths)
    print "BWA fail ", bwa_problems
    print "Warning: works only without GrowEnds()"
    [wrong_gaps, wrong_gap_empty, wrong_start, wrong_start_empty, wrong_end, wrong_end_empty, unknown] = cnt_wronglyclosedgaps(truepaths, alignedpaths, K) 
    print "Gaps stage problems (in, prefix, suffix, unknown) ", [wrong_gaps, wrong_gap_empty, wrong_start, wrong_start_empty, wrong_end, wrong_end_empty, unknown]
    [med_prefix, med_suffix, med_sum] = cnt_median_alignment_length(truepaths, alignedpaths, K)
    print "Median proportion of length of alignment between two fathest bwa hits to read length", med_prefix, med_suffix, med_sum 

    
    row_names = ["Total number of reads", \
                 "Read aligned to ref (#reads)",\
                 "Read aligned to ref and corresponding ref subseq to graph (#reads)",\
                 "Mapped with GAligner (#reads)",\
                 "Path is not equal to true path (#reads)",\
                 "Path is wrong. BWA hits uncertainty (#reads)",\
                 "Resulting BWA hits failure (#reads)", \
                 "Gap stage failure (#reads)", \
                 "Incorrect prefix/suffix (#reads)" , \
                 "Median length(in nucs) of skipped prefix/suffix/both"
                 ]
    res[fl] = {"Total number of reads" : len(reads), \
                 "Read aligned to ref (#reads)": len(reads) - cnt_badideal(reads, aligned_onref),\
                 "Read aligned to ref and corresponding ref subseq to graph (#reads)": len(reads) - cnt_badideal(reads, truepaths),\
                 "Mapped with GAligner (#reads)" : len(reads) - cnt_badideal(reads, truepaths) - cnt_notmapped(reads, truepaths, alignedpaths),\
                 "Path is not equal to true path (#reads)": path_problems,\
                 "Path is wrong. BWA hits uncertainty (#reads)": unknown - bwa_problems,\
                 "Resulting BWA hits failure (#reads)": bwa_problems, \
                 "Gap stage failure (#reads)": str(wrong_gaps) + " + " + str(wrong_gap_empty) + "(didn't closed)", \
                 "Incorrect prefix/suffix (#reads)" : str(wrong_start + wrong_start_empty) + "/" + str(wrong_end + wrong_end_empty), \
                 "Median length(in nucs) of skipped prefix/suffix/both" : \
                 str(med_prefix) + "/" + str(med_suffix) + "/" + str(med_sum)}



caption_below = ["Read aligned to ref and corresponding ref subseq to graph -- read aligned to reference by BWA MEM and alignment length > 0.8*(read length). After ref subseq mapped to graph -- the result of it is a true path.",\
                 "True path -- path produced by MapRead(), aligning sequence from reference that represents read",\
                 "Resulting BWA hits -- BWA hits after filtering, sorting etc., just before Gap closing stage", \
                 "Resulting BWA hits failure -- BWA hits are on edges that are not in true path or not in correct order",\
                 "Gap stage failure -- gap between two neibouring BWA hits wasn't closed by correct list of edges",\
                 "Prefix -- subpath of edges before first BWA hit edge in the path",\
                 "Suffix -- subpath of edges after last BWA hit edge in the path"]

table = make_table(res, row_names, caption_below, html_name[:-5])
save_html(table, "path_statistics_" + html_name)








