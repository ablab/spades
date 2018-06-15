import edlib

def edist(lst):
    if len(str(lst[0])) == 0:
        return len(str(lst[1]))
    if len(str(lst[1])) == 0:
        return len(str(lst[0]))
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
            for x in path_dirty.split("]")[:-1]:
                if x.startswith(","):
                    path.append(x.split()[1])
                else:
                    path.append(x.split()[0])

            res[cur_read] = { "len": rlen, "path": path, "edgelen": [int(x) for x in edgelen.split(",")[:-1]], "mapped": mapped_seq }
        fin.close()
        return res

    def load_alignments(self, filename):
        res = {}
        fin = open(filename, "r")
        bwa_num = 0
        cnt_empty = 0
        for ln in fin.readlines():
            cur_read, seq_starts, seq_ends, e_starts, e_ends, rlen, path, edgelen, bwa_path_dirty = ln.strip().split("\t")
            cur_read = cur_read.split(" ")[0]
            bwa_path = []
            edgelen_lst =[int(x) for x in edgelen.replace(";", "").split(",")[:-1]]
            empty = path.count(";")
            bwa_path_dirty = bwa_path_dirty.replace(";", "")
            path = path.replace(";", "")
            ranges = []
            edge_ranges = []
            for x in bwa_path_dirty.split("]")[:-1]:
                bwa_num += 1
                if x.startswith(","):
                    bwa_path.append(x.split()[1])
                else:
                    bwa_path.append(x.split()[0])
                ranges.append({"start": int(x.split("[")[1].split(",")[0]), "end": int(x.split("[")[1].split(",")[1])})
                edge_ranges.append({"start": int(x.split("(")[1].split(",")[0]), "end": int(x.split("(")[1].split(",")[1].split(")")[0])})
            initial_s = [int(x) for x in seq_starts.split(",")[:-1]] if "," in seq_starts else [int(seq_starts)]
            initial_e = [int(x) for x in seq_ends.split(",")[:-1]] if "," in seq_ends else [int(seq_ends)]
            mapped_s = [int(x) for x in e_starts.split(",")[:-1]] if "," in e_starts else [int(e_starts)]
            mapped_e = [int(x) for x in e_ends.split(",")[:-1]] if "," in e_ends else [int(e_ends)]
            mapped_len = 0
            for i in xrange(len(initial_s)):
                mapped_len += initial_e[i] - initial_s[i]
            res[cur_read] = { "len": rlen, "path": path.split(",")[:-1], "bwa_path": bwa_path, "edgelen": edgelen.split(",")[:-1], \
                                "initial_s": initial_s, "initial_e": initial_e, 
                                "mapped_s": mapped_s, "mapped_e": mapped_e, 
                                "mapped_len": mapped_len, 
                                "seq_ranges": ranges, "edge_ranges": edge_ranges, "empty": empty - 1}
            if empty - 1 > 0:
                cnt_empty += 1
        fin.close()
        print "Number of edges detected by bwa:", bwa_num
        print "gapped=", cnt_empty
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

    def divide_paths(self):
        alignedsubpath = {}
        badlyaligned = {}
        for r in self.alignedpaths.keys():
            if r not in self.truepaths.keys():
                continue
            if self.truepaths[r]["path"] != self.alignedpaths[r]["path"] or self.alignedpaths[r]["empty"] > 0:
                if ",".join(self.alignedpaths[r]["path"]) in ",".join(self.truepaths[r]["path"]) and self.alignedpaths[r]["empty"] == 0:
                    alignedsubpath[r] = self.alignedpaths[r]
                elif self.truepaths[r]["path"] != self.alignedpaths[r]["path"] or self.alignedpaths[r]["empty"] > 0:
                    badlyaligned[r] = self.alignedpaths[r]
                if ",".join(self.alignedpaths[r]["path"]) in ",".join(self.truepaths[r]["path"]) and self.alignedpaths[r]["empty"] > 0:
                    print "Ranges_failure readname=",r 
        return [alignedsubpath, badlyaligned]


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
            for e in self.alignedpaths[r]["bwa_path"]:
                if e in self.truepaths[r]["path"]:
                    j += 1
                else:
                    bwa_fails.append(e)
            if j < len(self.alignedpaths[r]["bwa_path"]) :
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
        path = self.alignedpaths[r]["path"]
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
                print ",".join(self.alignedpaths[r]["bwa_path"])
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
                print ",".join(self.alignedpaths[r]["bwa_path"])
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
        ln = 0
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
        ln = 0
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

        med_prefix = 0
        if len(res_prefix) > 0:
            med_prefix = sorted(res_prefix)[len(res_prefix)/2]

        med_suffix = 0
        if len(res_suffix) > 0:
            med_suffix = sorted(res_suffix)[len(res_suffix)/2]

        med_plus = 0
        if len(res_plus) > 0:
            med_plus = sorted(res_plus)[len(res_plus)/2]
        return {"prefix_len": med_prefix, "suffix_len": med_suffix, "sum_len": med_plus }


class GapEditDistanceCounter:

    def __init__(self, bwahitsMapper, edges, K):
        self.bwahits_mapper = bwahitsMapper
        self.K = K
        self.reads = bwahitsMapper.reads
        self.truepaths = bwahitsMapper.truepaths
        self.alignedpaths = bwahitsMapper.alignedpaths
        self.edges = edges


    def extract_subpath(self, subpath, range1, range2):
        res = self.edges[subpath[0]][range1["end"]: len(self.edges[subpath[0]]) - self.K]
        if range1["end"] > len(self.edges[subpath[0]]) - self.K:
            print "WARNING"
        for ind in xrange(1, len(subpath) - 1):
            res += self.edges[subpath[ind]][: len(self.edges[subpath[ind]]) - self.K]
        res += self.edges[subpath[len(subpath) - 1]][ : range2["start"]]
        return res


    def extract_subseq(self, r, range1, range2):
        seq = self.reads[r]
        return seq[range1["end"]: range2["start"]]


    def count_gap_ed(self):
        cnt = 0
        total = 0
        for r in self.alignedpaths.keys():
            if r not in self.truepaths:
                continue
            if self.bwahits_mapper.is_unique(r):
                true_ind = self.bwahits_mapper.get_bwa_inds_true(r)
                aligned_ind = self.bwahits_mapper.get_bwa_inds_aligned(r)
                for j in xrange(2, len(true_ind) - 1):
                    if self.truepaths[r]["path"][true_ind[j - 1]: true_ind[j]] != self.alignedpaths[r]["path"][aligned_ind[j - 1]: aligned_ind[j]]:
                        galigner_seq = self.extract_subpath(self.alignedpaths[r]["path"][aligned_ind[j - 1]: aligned_ind[j] + 1],\
                                                    self.alignedpaths[r]["edge_ranges"][j - 2], self.alignedpaths[r]["edge_ranges"][j - 1])
                        true_seq = self.extract_subpath(self.truepaths[r]["path"][true_ind[j - 1]: true_ind[j] + 1],\
                                                    self.alignedpaths[r]["edge_ranges"][j - 2], self.alignedpaths[r]["edge_ranges"][j - 1])
                        read_seq = self.extract_subseq(r, self.alignedpaths[r]["seq_ranges"][j - 2], self.alignedpaths[r]["seq_ranges"][j - 1])
                        if edist([read_seq, true_seq]) < edist([read_seq, galigner_seq]):
                            cnt += 1
                        total += 1
                        # if edist([read_seq, true_seq]) < edist([read_seq, galigner_seq]) and len(self.alignedpaths[r]["path"][aligned_ind[j - 1]: aligned_ind[j]+1]) > 2: 
                        #     print r, edist([read_seq, true_seq]), edist([read_seq, galigner_seq]), edist([galigner_seq, true_seq]),len(read_seq)
                        #     print self.truepaths[r]["path"][true_ind[j - 1]: true_ind[j] + 1]
                        #     print self.alignedpaths[r]["path"][aligned_ind[j - 1]: aligned_ind[j] + 1]
                        #     print ""
        return {"failed_ed": cnt, "ed": total} 
