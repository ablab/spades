import sys
import re

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

logname = sys.argv[1]
idealfile = sys.argv[2]
prefix = sys.argv[3]

def load_la_paths(logname, prefix):
    edge_sets = {}
    with open(logname, "r") as fin:
        s = fin.readline()
        pattern = re.compile("ReadName=S\d_\d+ edge=\d+ edge_len=\d+ e_start=\d+ e_end=\d+")
        while s != "":
            res = re.findall(pattern, s.strip())
            for it in res:
                readname = it.split(" ")[0][len("ReadName="):]
                edgeid = it.split(" ")[1][len("edge="):]
                edge_len = int(it.split(" ")[2][len("edge_len="):])
                e_start = int(it.split(" ")[3][len("e_start="):])
                e_end = int(it.split(" ")[4][len("e_end="):])
                if readname not in edge_sets.keys():
                    edge_sets[readname] = []
                edge_sets[readname].append([edgeid, str(e_end - e_start), e_start*1.0/edge_len > 0.3 and (edge_len - e_end) > 0.3])
            s = fin.readline()

    with open(prefix+"_la.tsv", "w") as fout:
        for read in sorted(edge_sets.keys()):
            fout.write(read + "\t" + ",".join([y[0] + "(" + y[1] + "/" + str(y[2]) + ")" for y in sorted(edge_sets[read], key=lambda x: x[0])]) + "\n")

    return edge_sets

def load_ideal_paths(idealfile, prefix):
    edge_sets = {}
    with open(idealfile, "r") as fin:
        for ln in fin.readlines():
            [name, nuc_len, path_len, path, edge_len, aln] = ln.strip().split("\t")
            if name not in edge_sets:
                edge_sets[name] = []
            path_lst = [x.split(" (")[0] for x in path.split(") ")[:-1]]
            edge_len_lst = edge_len.split(",")
            for i in xrange(len(path_lst)):
                edge_sets[name].append([path_lst[i], edge_len_lst[i]])

    with open(prefix + "_ideal.tsv", "w") as fout:
        for read in sorted(edge_sets.keys()):
            fout.write(read + "\t" + ",".join([y[0] + "(" + y[1] + ")" for y in sorted(edge_sets[read], key=lambda x: x[0]) ]) + "\n")

    return edge_sets


ideal_paths = load_ideal_paths(idealfile, prefix)
la_paths = load_la_paths(logname, prefix)
cnt = 0
cnt_in = 0
only_ideal_lens = []
only_la_lens = []
total_alone = 0
fp_alone = 0
with open(prefix + "_dif.tsv", "w") as fout:
    for k in sorted(la_paths.keys()):
        ideal_set = set([x[0] for x in ideal_paths[k]])
        la_set = set([x[0] for x in la_paths[k]])
        la_alone_set = set()
        for x in la_paths[k]:
            if x[2]:
                la_alone_set.add(x[0])
        total_alone += len(la_alone_set)
        if ideal_set != la_set:
            s = ""
            if ideal_set.issubset(la_set):
                cnt_in += 1
                s = "IN"
            else:
                s = "OUT"
            only_ideal = ideal_set - la_set
            for it in ideal_paths[k]:
                if it[0] in only_ideal:
                    only_ideal_lens.append(int(it[1]))

            alone_la = la_alone_set - ideal_set
            fp_alone += len(alone_la)
            only_la = la_set - ideal_set
            for it in la_paths[k]:
                if it[0] in only_la:
                    only_la_lens.append(int(it[1]))

            fout.write(k + "\t" + ",".join([y[0] + "(" + y[1] + ")" for y in sorted(ideal_paths[k], key=lambda x: x[0])])
                       + "\t" + ",".join([y[0] + "(" + y[1] + ")" for y in sorted(la_paths[k], key=lambda x: x[0])]) + "\t" + s + "\n")
            cnt += 1

print cnt_in, cnt, len(ideal_paths.keys()), len(la_paths.keys())
print len(only_ideal_lens), len(only_la_lens)
print "FN ", sorted(only_ideal_lens)[len(only_ideal_lens)/2], sum(only_ideal_lens)*1.0/len(only_ideal_lens), sorted(only_ideal_lens)[len(only_ideal_lens) - 1]
print "FP ", sorted(only_la_lens)[len(only_la_lens)/2], sum(only_la_lens)*1.0/len(only_la_lens), sorted(only_la_lens)[len(only_la_lens) - 1]
print "FP alone ", fp_alone*1.0/total_alone
print "Draw FN.."
plt.figure()
plt.hist(only_ideal_lens, bins = 50)
plt.title('FN')
plt.xlabel('mapping length')
#plt.ylabel('Another Label')
plt.savefig(prefix + "_hist_ideal.png")
print "Draw FP.."
plt.figure()
plt.hist(only_la_lens, bins=50)
plt.title('FP')
plt.xlabel('mapping length')
#plt.ylabel('Another Label')
plt.savefig(prefix + "_hist_la.png")
