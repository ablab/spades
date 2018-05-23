import sys
import yaml
import edlib

def edist(lst):
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

class DataLoader:

    def load(self, filename, datatype):
        if datatype == "fasta":
            return self.load_reads(filename)
        elif datatype == "true_paths":
            return self.load_truepaths(filename)
        elif datatype == "galigner_paths":
            return self.load_alignments(filename)
        elif datatype == "txt":
            return self.load_list(filename)

    def load_list(self, filename):
        res = []
        with open(filename, "r") as fin:
            res = [ln.strip() for ln in fin.readlines()]

        res = sorted(res, key= lambda x: x[::-1])
        return res

    def load_reads(self, filename):
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

    def load_truepaths(self, filename):
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

    def load_alignments(self, filename):
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
                                "mapped_s":int(seq_start), "mapped_e":int(seq_end), "seq": seq , "seq_end": seq_ends, "empty": empty - 1}
        fin.close()
        print "Number of edges detected by bwa:", bwa_num
        return res


class GeneralStatisticsCounter:

    def __init__(self, reads, truepaths, alignedpaths):
        self.reads = reads
        self.truepaths = truepaths
        self.alignedpaths = alignedpaths

    def cnt_badideal(self):
        res = 0
        for r in self.reads.keys():
            if r not in self.truepaths:
                res += 1
        return res

    def cnt_notmapped(self):
        res = 0
        for r in self.reads.keys():
            if r not in self.truepaths:
                continue
            if r not in self.alignedpaths:
                #print r
                res += 1
        return res

    def cnt_problempaths(self):
        res = 0
        for r in self.alignedpaths.keys():
            if r not in self.truepaths.keys():
                continue
            if self.truepaths[r]["path"] != self.alignedpaths[r]["path"] or self.alignedpaths[r]["empty"] > 0:
                res += 1
        return res

class BWAhitsMapper:

    def __init__(self, reads, truepaths, alignedpaths):
        self.reads = reads
        self.truepaths = truepaths
        self.alignedpaths = alignedpaths

    def cnt_problembwa(self, print_fails = False):
        res = 0
        total_path = 0
        for r in self.alignedpaths.keys():
            if r not in self.truepaths:
                continue
            j = 0
            total_path += 1
            # for e in truepaths[r]["path"]:
            #     if j < len(alignedpaths[r]["bwa_path"]) and e == alignedpaths[r]["bwa_path"][j]:
            #         j += 1
            bwa_fails = []
            for e in alignedpaths[r]["bwa_path"]:
                if e in truepaths[r]["path"]:
                    j += 1
                else:
                    bwa_fails.append(e)
            if j < len(alignedpaths[r]["bwa_path"]) :
                res += 1
                if print_fails:
                   print "BWA_Fail readname=", r, " failed hits:", ",".join([str(x) for x in bwa_fails])
        return res

    def get_bwa_inds_true(self, r):
        path = self.truepaths[r]["path"]
        bwapath = self.alignedpaths[r]["bwa_path"]
        ind = []
        j = 0
        for i in xrange(len(path)):
            if j < len(bwapath) and path[i] == bwapath[j]:
                j += 1
                ind.append(i)
        ind.insert(0, 0)
        ind.append(len(path))
        return ind

    def get_bwa_inds_aligned(self, r):
        path = self.alignedpaths[r]["edge_tag"]
        ind = []
        for i in xrange(len(path)):
            if path[i] == "bwa":
                ind.append(i)
        ind.insert(0, 0)
        ind.append(len(path))
        return ind

    def is_unique(self, r):
        path = self.truepaths[r]["path"]
        subpath = self.alignedpaths[r]["bwa_path"]
        ind = self.get_bwa_inds_true(r)
        if len(ind) < len(subpath) + 2:
            return False
        for j in xrange(1, len(subpath) + 1):
            if subpath[j - 1] in path[ind[j - 1]: ind[j]] or subpath[j - 1] in path[ind[j] + 1: ind[j + 1]]:
                return False
        return True


class GapsStatistics:

    def __init__(self, bwahitsMapper, edges, K, print_fails = False):
        self.bwahitsMapper = bwahitsMapper
        self.edges = edges
        self.K = K
        self.print_fails = print_fails
        self.reads = bwahitsMapper.reads
        self.truepaths = bwahitsMapper.truepaths
        self.alignedpaths = bwahitsMapper.alignedpaths

    def is_wrong_start(self, r, true_ind, ind):
        empty = False
        truepath = self.truepaths[r]["path"]
        path = self.alignedpaths[r]["path"]
        edgelen = self.truepaths[r]["edgelen"]
        aedgelen = self.alignedpaths[r]["edgelen"]
        if truepath[true_ind[0]: true_ind[1]] == path[ind[0]: ind[1]] :
            return [False, empty]
        else:
            if self.print_fails:
                print "Wrong_start readname=", r
                print ",".join(truepath[true_ind[0]: true_ind[1]]), ",".join([str(x) for x in edgelen[true_ind[0]: true_ind[1]]])
                print ",".join(path[ind[0]: ind[1]]), ",".join(aedgelen[ind[0]: ind[1]])
                print ""
            i = ind[1]
            j = true_ind[1]
            while i >=0 and j >=0 and path[i] == truepath[j]:
                i -= 1
                j -= 1
            if len(edgelen[:j + 1]) == 0:
                if sum([int(x) for x in aedgelen[:i + 1]]) <= self.K:    
                    return [False, empty]
            if len(path[ind[0]: ind[1]]) == 0:
                empty = True
            return [True, empty]

    def is_wrong_end(self, r, true_ind, ind):
        empty = False
        truepath = self.truepaths[r]["path"]
        path = self.alignedpaths[r]["path"]
        edgelen = self.truepaths[r]["edgelen"]
        aedgelen = self.alignedpaths[r]["edgelen"]
        if truepath[true_ind[-2]: true_ind[-1]] == path[ind[-2]: ind[-1]] :
            return [False, empty]
        else:
            if self.print_fails:
                print "Wrong_end readname=", r
                print ",".join(truepath[true_ind[-2]: true_ind[-1]]), ",".join([str(x) for x in edgelen[true_ind[-2]: true_ind[-1]]])
                print ",".join(path[ind[-2]: ind[-1]]), ",".join(aedgelen[ind[-2]: ind[-1]])
                print ""
            i = ind[-2]
            j = true_ind[-2]
            while i < len(path) and j < len(truepath) and path[i] == truepath[j]:
                i += 1
                j += 1
            if len(edgelen[j:]) == 0:
                if sum([int(x) for x in aedgelen[i:]]) <= self.K:  
                    return [False, empty]
            if len(path[ind[-2]: ind[-1]]) == 1:
                empty = True
            return [True, empty]

    def restore_path(self, lst):
        return "".join([self.edges[it][:-self.K] for it in lst])

    def is_wrong_gap(self, r, true_ind, ind):
        wrong_gap = 0
        truepath = self.truepaths[r]["path"]
        path = self.alignedpaths[r]["path"]
        empty = self.alignedpaths[r]["empty"]
        true_lens = self.truepaths[r]["edgelen"]
        lens = self.alignedpaths[r]["edgelen"]
        for j in xrange(2, len(true_ind) - 1):
            if truepath[true_ind[j - 1]: true_ind[j]] != path[ind[j - 1]: ind[j]]:
                wrong_gap += 1
                if self.print_fails:
                    if len(path[ind[j - 1]: ind[j]]) > 1:
                        print r
                        print ",".join(truepath[true_ind[j - 1]+1: true_ind[j]]), sum([int(x) for x in true_lens[true_ind[j - 1]+1: true_ind[j]]])
                        len1 = sum([int(x) for x in lens[ind[j - 1]+1: ind[j]]])
                        len2 = sum([int(x) for x in true_lens[true_ind[j - 1]+1: true_ind[j]]])
                        print ",".join(path[ind[j - 1]+1: ind[j]]), sum([int(x) for x in lens[ind[j - 1]+1: ind[j]]])
                        print "ed=", edist([self.restore_path(truepath[true_ind[j - 1] + 1: true_ind[j]]), self.restore_path(path[ind[j - 1] + 1: ind[j]]) ]), len(self.restore_path(truepath[true_ind[j - 1] + 1: true_ind[j]]))
                        if edist([self.restore_path(truepath[true_ind[j - 1] + 1: true_ind[j]]), self.restore_path(path[ind[j - 1] + 1: ind[j]]) ]) >= 10:
                            print "Big",abs(len1 - len2), abs(len1 - len2)*100/len2
                        print ""
                    else:
                        print "Empty: ", r, " len=", sum([int(x) for x in true_lens[true_ind[j - 1]+1: true_ind[j]]])
                        print ""
        return wrong_gap - empty > 0


    def cnt_wronglyclosedgaps(self):
        truepaths = self.truepaths
        alignedpaths = self.alignedpaths
        res, res_start, res_end = 0, 0, 0
        res_empty, res_start_empty, res_end_empty = 0, 0, 0
        unknown = 0
        for r in alignedpaths.keys():
            if r not in truepaths:
                continue
            if not self.bwahitsMapper.is_unique(r):
                unknown += 1
            else:
                true_ind = self.bwahitsMapper.get_bwa_inds_true(r)
                aligned_ind = self.bwahitsMapper.get_bwa_inds_aligned(r)
                has_wrong_start, empty = self.is_wrong_start(r, true_ind, aligned_ind)
                if has_wrong_start:
                    if empty:
                        res_start_empty += 1
                    else:
                        res_start += 1
                has_wrong_end, empty = self.is_wrong_end(r, true_ind, aligned_ind)
                if has_wrong_end:
                    if empty:
                        res_end_empty += 1
                    else:
                        res_end += 1
                has_wrong_gap = self.is_wrong_gap(r, true_ind, aligned_ind)
                if has_wrong_gap:
                    res += 1
                if alignedpaths[r]["empty"] > 0:
                    res_empty += 1
        return {"wrong_filled_gaps": res, "empty_gaps": res_empty, "wrong_filled_start": res_start, "empty_start": res_start_empty, "wrong_filled_end": res_end, "empty_end": res_end_empty, "unknown_num": unknown}

class GapsLengthStatistics:

    def __init__(self, bwahitsMapper, K):
        self.bwahitsMapper = bwahitsMapper
        self.K = K
        self.reads = bwahitsMapper.reads
        self.truepaths = bwahitsMapper.truepaths
        self.alignedpaths = bwahitsMapper.alignedpaths

    def is_wrong_start(self, r, true_ind, ind):
        empty = False
        truepath = self.truepaths[r]["path"]
        path = self.alignedpaths[r]["path"]
        edgelen = self.truepaths[r]["edgelen"]
        aedgelen = self.alignedpaths[r]["edgelen"]
        ln = -1
        if truepath[true_ind[0]: true_ind[1]] == path[ind[0]: ind[1]] :
            return [False, empty, ln]
        else:
            i = ind[1]
            j = true_ind[1]
            while i >=0 and j >=0 and path[i] == truepath[j]:
                i -= 1
                j -= 1
            if len(edgelen[:j + 1]) == 0:
                if sum([int(x) for x in aedgelen[:i + 1]]) <= self.K:    
                    return [False, empty, ln]
            else:
                ln = 0
                for it in edgelen[:j + 1]:
                    ln += it
            if len(path[ind[0]: ind[1]]) == 0:
                empty = True
            return [True, empty, ln]

    def is_wrong_end(self, r, true_ind, ind):
        empty = False
        truepath = self.truepaths[r]["path"]
        path = self.alignedpaths[r]["path"]
        edgelen = self.truepaths[r]["edgelen"]
        aedgelen = self.alignedpaths[r]["edgelen"]
        ln = -1
        if truepath[true_ind[-2]: true_ind[-1]] == path[ind[-2]: ind[-1]] :
            return [False, empty, ln]
        else:
            i = ind[-2]
            j = true_ind[-2]
            while i < len(path) and j < len(truepath) and path[i] == truepath[j]:
                i += 1
                j += 1
            if len(edgelen[j:]) == 0:
                if sum([int(x) for x in aedgelen[i:]]) <= self.K:  
                    return [False, empty, ln]
            else:
                ln = 0
                for it in edgelen[j:]:
                    ln += it
            if len(path[ind[-2]: ind[-1]]) == 1:
                empty = True
            return [True, empty, ln]


    def cnt_median_alignment_length(self):
        unknown = 0
        tpaths = self.truepaths
        apaths = self.alignedpaths
        res_prefix = []
        res_suffix = []
        res_plus = []
        wrong_start = 0
        wrong_end = 0
        for r in apaths.keys():
            if r not in tpaths:
                continue
            if not self.bwahitsMapper.is_unique(r):
                unknown += 1
            else:
                true_ind = self.bwahitsMapper.get_bwa_inds_true(r)
                aligned_ind = self.bwahitsMapper.get_bwa_inds_aligned(r)
                has_wrong_start, empty, ln = self.is_wrong_start(r, true_ind, aligned_ind)
                m = 0
                if has_wrong_start:
                    wrong_start += 1
                    res_prefix.append(ln)
                    m += ln
                has_wrong_end, empty, ln = self.is_wrong_end(r, true_ind, aligned_ind)
                if has_wrong_end:
                    wrong_end += 1
                    res_suffix.append(ln)
                    m += ln
    
                if has_wrong_end or has_wrong_start:
                    res_plus.append(m)

        return {"prefix_len": sorted(res_prefix)[len(res_prefix)/2], "suffix_len": sorted(res_suffix)[len(res_suffix)/2], "sum_len": sorted(res_plus)[len(res_plus)/2] }

def make_table(results, row_names, caption, name):
    html = """<html><table border="1"><caption>{}</caption><tr><th></th>""".format(name)
    for run_name in sorted(results.keys(), key= lambda x: x[::-1]):
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

def get_name(path):
    name = path.split("/")[-1]
    res = ""
    if name.startswith("dima_filtering"):
        res = "branch: new_weights; "
    else:
        res = "branch: master; "
    if "_bf_" in name:
        res += "brute_force; "
    else:
        res += "dijkstra; "
    if "_ends" in name:
        res += "restore_ends; "
    if "_ideal_" in name:
        res += "ideal_reads; "

    return res

def save_html(s, fl):
    with open(fl, "w") as fout:
        fout.write(s)


with open(sys.argv[1], 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

data_loader = DataLoader()
reads = data_loader.load(cfg["reads_path"], "fasta")
edges = data_loader.load(cfg["edges"], "fasta")
truepaths = data_loader.load(cfg["truepaths_tsv"], "true_paths")
aligned_files = data_loader.load(cfg["galigner_out"], "txt")
print "\n".join(aligned_files)
K = int(cfg["k"])
html_name = cfg["print_html"]

res = {}
for fl in aligned_files:
    alignedpaths = data_loader.load(fl, "galigner_paths")

    general_stats = GeneralStatisticsCounter(reads, truepaths, alignedpaths)
    badideal = general_stats.cnt_badideal()
    notmapped = general_stats.cnt_notmapped() 
    path_problems = general_stats.cnt_problempaths()

    bwahits_mapper = BWAhitsMapper(reads, truepaths, alignedpaths)
    bwa_problems = bwahits_mapper.cnt_problembwa(cfg["print_bwa_hits_failure"])

    gaps_statistics = GapsStatistics(bwahits_mapper, edges, K, cfg["print_gaps_failure"])
    gaps_cnt_stats = gaps_statistics.cnt_wronglyclosedgaps() 
    
    gaps_stat_len = GapsLengthStatistics(bwahits_mapper, K)
    gaps_len_stats = gaps_stat_len.cnt_median_alignment_length()
     
    print "Total=", len(reads), " ideal=", len(reads) - badideal,  " notmapped=", notmapped
    print "Paths with problems ", path_problems
    print "BWA fail ", bwa_problems
    print "Gaps stage problems (in, prefix, suffix, unknown) ", gaps_cnt_stats
    print "Median proportion of length of alignment between two fathest bwa hits to read length", gaps_len_stats

    def make_str(cnt, total):
        return str(cnt) + " (" + str(cnt*100/total) + "%)"

    def make_str2(cnt, cnt_empty, total):
        return str(cnt) + " + " + str(cnt_empty) + " (" + str((cnt + cnt_empty)*100/total) + "%)"

    res[fl.split("/")[-1]] = {
                 "Total number of reads": str(len(reads) - badideal)+ " (100%)",\
                 "Mapped with GAligner (#reads)" : make_str(len(reads) - badideal - notmapped, len(reads) - badideal),\
                 "Path is not equal to true path (#reads)": make_str(path_problems, len(reads) - badideal),\
                 "Path is wrong. BWA hits uncertainty (#reads)": make_str(gaps_cnt_stats["unknown_num"] - bwa_problems, len(reads) - badideal) ,\
                 "Resulting BWA hits failure (#reads)": make_str(bwa_problems, len(reads) - badideal), \
                 "Gap stage failure (#reads)": make_str2(gaps_cnt_stats["wrong_filled_gaps"], gaps_cnt_stats["empty_gaps"], len(reads) - badideal), \
                 "Incorrect prefix/suffix (#reads)" : make_str2(gaps_cnt_stats["wrong_filled_start"], gaps_cnt_stats["empty_start"], len(reads) - badideal) + " / " +\
                                                      make_str2(gaps_cnt_stats["wrong_filled_end"], gaps_cnt_stats["empty_end"], len(reads) - badideal), \
                 "Median length(in nucs) of skipped prefix/suffix/both" : \
                 str(gaps_len_stats["prefix_len"]) + "/" + str(gaps_len_stats["suffix_len"]) + "/" + str(gaps_len_stats["sum_len"])}

row_names = [
                 "Total number of reads",\
                 "Mapped with GAligner (#reads)",\
                 "Path is not equal to true path (#reads)",\
                 "Path is wrong. BWA hits uncertainty (#reads)",\
                 "Resulting BWA hits failure (#reads)", \
                 "Gap stage failure (#reads)", \
                 "Incorrect prefix/suffix (#reads)" , \
                 "Median length(in nucs) of skipped prefix/suffix/both"
                 ]

caption_below = ["Read aligned to ref and corresponding ref subseq to graph -- read aligned to reference by BWA MEM and alignment length > 0.8*(read length). After ref subseq mapped to graph -- the result of it is a true path.",\
                 "True path -- path produced by MapRead(), aligning sequence from reference that represents read",\
                 "Resulting BWA hits -- BWA hits after filtering, sorting etc., just before Gap closing stage", \
                 "Resulting BWA hits failure -- BWA hits are on edges that are not in true path or not in correct order",\
                 "Gap stage failure -- gap between two neibouring BWA hits wasn't closed by correct list of edges",\
                 "Prefix -- subpath of edges before first BWA hit edge in the path",\
                 "Suffix -- subpath of edges after last BWA hit edge in the path"]

table = make_table(res, row_names, caption_below, html_name[:-5])
save_html(table, "path_statistics_" + html_name)








