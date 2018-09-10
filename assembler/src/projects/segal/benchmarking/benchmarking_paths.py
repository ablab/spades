import edlib
import pandas as pd
import matplotlib
import sys
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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
        res = []
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

            res.append({"r_name": cur_read, "len": rlen, "true_path": ",".join(path), "truemapped_seq": mapped_seq })
        fin.close()
        return res

    def load_segal_paths(self, filename, reads):
        res = []
        fin = open(filename, "r")
        for ln in fin.readlines():
            cur_read, seq_starts, seq_ends, e_starts, e_ends, rlen, path, edgelen, bwa_path_dirty, seqs = ln.strip().split("\t")
            cur_read = cur_read.split(" ")[0]
            initial_s = [int(x) for x in seq_starts.split(",")[:-1]] if "," in seq_starts else [int(seq_starts)]
            initial_e = [int(x) for x in seq_ends.split(",")[:-1]] if "," in seq_ends else [int(seq_ends)]
            mapped_s = [int(x) for x in e_starts.split(",")[:-1]] if "," in e_starts else [int(e_starts)]
            mapped_e = [int(x) for x in e_ends.split(",")[:-1]] if "," in e_ends else [int(e_ends)]
            max_ind = 0
            for i in xrange(1, len(initial_s)):
                if initial_e[max_ind] - initial_s[max_ind] < initial_e[i] - initial_s[i]:
                    max_ind = i
            if initial_e[max_ind] - initial_s[max_ind] + 1 > 0.*len(reads[cur_read]):
                res.append({"r_name": cur_read, "mapping_len": initial_e[max_ind] - initial_s[max_ind] + 1, \
                                "path": path.split(";")[max_ind][:-1],
                                "s_start": initial_s[max_ind], "s_end": initial_e[max_ind], \
                                "mapped_seq": seqs.split(";")[max_ind]})
        fin.close()
        return res


def check_ed(seq, x, tpath, gapath):
    tmapped = x["truemapped_seq"]
    gamapped = x["mapped_seq"]
    s, e = x["s_start"], x["s_end"]
    # if tpath == gapath:
    #     if s != 0 or e != len(seq):
    #         print x["r_name"]
    #     elif edist([seq, tmapped]) < edist([seq, gamapped]):
    #         print "More:", x["r_name"], edist([seq, tmapped]), edist([seq, gamapped])

    if s != 0 or e != len(seq):
        return False
    if edist([seq, tmapped]) > edist([seq, gamapped]):
        return False
    return True 

def check_path(seq, x):
    tpath = x["true_path"]
    gapath = x["path"]
    s, e = x["s_start"], x["s_end"]
    if tpath == gapath: # and s == 0 and e == len(seq):
        return True
    return False


def print_stats(reads, ppaths, gapaths):
    perfect_df = pd.DataFrame(ppaths)
    ga_df = pd.DataFrame(gapaths)
    df = pd.merge(perfect_df, ga_df, on="r_name")
    df["better_ed"] = df.apply(lambda x: check_ed(reads[x["r_name"]], x, x["true_path"], x["path"]), axis = 1)
    df["equal_paths"] = df.apply(lambda x: check_path(reads[x["r_name"]], x), axis = 1)
    tp_baded = df.apply(lambda x: (not x["better_ed"]) and x["equal_paths"], axis = 1).sum()
    tp_gooded = df.apply(lambda x: x["better_ed"] and x["equal_paths"], axis = 1).sum()
    print "TruePath:", "BadED=", tp_baded, \
                      "GoodED=", tp_gooded
    print "WrongPath:", "BadED=", df.apply(lambda x: (not x["better_ed"]) and (not x["equal_paths"]), axis = 1).sum(), \
                      "GoodED=", df.apply(lambda x: x["better_ed"] and (not x["equal_paths"]), axis = 1).sum()


    print (tp_baded + tp_gooded)*100/len(df), "%"
    print ""



if __name__ == "__main__":
    datapath = "/Sid/tdvorkina/gralign/benchmarking_paths"
    for org in ["ecoli", "celegans"]:
        print org
        for read_type in ["simnp", "simpb", "realpb", "realnp"]:
            print read_type
            org_path = datapath + "/" + org + "/"
            reads_file = org_path + "/input/" + read_type + "_mapped_10000.fasta"
            galigner_res_file = org_path + "/SeGal/output/aln_" + read_type + ".tsv"
            perfectpath_res_file = org_path + "/input/" + read_type + "_true_mapping.tsv"
            dl = DataLoader()
            reads = dl.load_reads(reads_file)
            galigner_res = dl.load_segal_paths(galigner_res_file, reads)
            perfectpaths = dl.load_truepaths(perfectpath_res_file)
            print_stats(reads, perfectpaths, galigner_res)
    # print "Draw reads len distribution.."
    # lens = [len(reads[k]) for k in reads.keys()
    # plt.hist(lens)
    # plt.xlabel("Read length")
    # plt.savefig("./len_dist.png")