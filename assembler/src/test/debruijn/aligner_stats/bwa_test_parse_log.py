import sys
import re

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

logname = sys.argv[1]
idealfile = sys.argv[2]
readsfile = sys.argv[3]
prefix = sys.argv[4]
thresholds = [0, 50, 100, 200, 300, 400, 500, 700, 900, 1000, 1500]

def load_reads_len(filename):
    res = {}
    key = ""
    cur_str = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip().split(" ")[0]
                cur_str = ""
                res[key] = 0
            else:
                cur_str += ln.strip()
                res[key] = len(cur_str)

    return res

def load_la_paths(logname, prefix, threshold, reads_len):
    edge_sets = {}
    with open(logname, "r") as fin:
        s = fin.readline()
        pattern = re.compile("ReadName=S\d_\d+ edge=\d+ edge_len=\d+ e_start=\d+ e_end=\d+ r_start=\d+ r_end=\d+")
        while s != "":
            res = re.findall(pattern, s.strip())
            for it in res:
                lst = it.split(" ")
                readname = lst[0][len("ReadName="):]
                edgeid = lst[1][len("edge="):]
                edge_len = int(lst[2][len("edge_len="):])
                e_start = int(lst[3][len("e_start="):])
                e_end = int(lst[4][len("e_end="):])
                r_start = int(lst[5][len("r_start="):])
                r_end = int(lst[6][len("r_end="):])
                if readname not in edge_sets and edge_len > threshold:
                    edge_sets[readname] = []
                if edge_len > threshold:
                    edge_sets[readname].append([edgeid, str(e_end - e_start), 2*(e_end - e_start) < min(r_start, e_start) + min(edge_len - e_end, reads_len[readname] - r_end)])
                #print [edgeid, str(e_end - e_start), ]
            s = fin.readline()

    with open(prefix+"_la.tsv", "w") as fout:
        for read in sorted(edge_sets.keys()):
            fout.write(read + "\t" + ",".join([y[0] + "(" + y[1] + "/" + str(y[2]) + ")" for y in sorted(edge_sets[read], key=lambda x: x[0])]) + "\n")

    return edge_sets

def load_ideal_paths(idealfile, prefix, threshold):
    edge_sets = {}
    with open(idealfile, "r") as fin:
        for ln in fin.readlines():
            [name, nuc_len, path_len, path, edge_len, whole_edge_len, aln] = ln.strip().split("\t")
            path_lst = [x.split(" (")[0] for x in path.split(") ")[:-1]]
            edge_len_lst = edge_len.split(",")[:-1]
            whole_edge_len_lst = [int(x) for x in whole_edge_len.split(",")[:-1]]
            add = False
            for e in whole_edge_len_lst:
                if e > threshold:
                    add = True
                    break
            if name not in edge_sets and add:
                edge_sets[name] = []
            if add:
                for i in xrange(len(path_lst)):
                    if whole_edge_len_lst[i] > threshold:
                        edge_sets[name].append([path_lst[i], edge_len_lst[i]])

    with open(prefix + "_ideal.tsv", "w") as fout:
        for read in sorted(edge_sets.keys()):
            fout.write(read + "\t" + ",".join([y[0] + "(" + y[1] + ")" for y in sorted(edge_sets[read], key=lambda x: x[0]) ]) + "\n")

    return edge_sets

precisions = []
recalls = []
f1scores = []
reads_len = load_reads_len(readsfile)

for threshold in thresholds:
    ideal_paths = load_ideal_paths(idealfile, prefix, threshold)
    la_paths = load_la_paths(logname, prefix, threshold, reads_len)

    cnt = 0
    cnt_in = 0
    only_ideal_lens = []
    only_la_lens = []
    total_alone = 0
    fp_alone = 0
    ideal_name_set = set(ideal_paths.keys())
    la_name_set = set(la_paths.keys())
    in_name_set = ideal_name_set.intersection(la_name_set)
    print len(ideal_name_set - la_name_set), len(ideal_name_set)
    right = 0
    right_small = 0
    right_0 = 0
    with open(prefix + "_" + str(threshold) + "_dif.tsv", "w") as fout:
        for k in sorted(ideal_name_set):
            ideal_set = set()
            if k in ideal_paths:
                ideal_set = set([x[0] for x in ideal_paths[k]])
            la_set = set()
            if k in la_paths:
                la_set = set([x[0] for x in la_paths[k]])
            la_alone_set = set()
            if k in la_paths:
                for x in la_paths[k]:
                    if x[2]:
                        la_alone_set.add(x[0])
                if len(ideal_set) > 1:
                    total_alone += len(la_alone_set)
            right += len(la_set.intersection(ideal_set))
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
                if len(ideal_set) > 1:
                    alone_la = la_alone_set - ideal_set
                    fp_alone += len(alone_la)
                only_la = la_set - ideal_set
                if k in la_paths:
                    for it in la_paths[k]:
                        if it[0] in only_la:
                            only_la_lens.append(int(it[1]))
                if k in la_paths:
                    fout.write(k + "\t" + ",".join([y[0] + "(" + y[1] + ")" for y in sorted(ideal_paths[k], key=lambda x: x[0])])
                               + "\t" + ",".join([y[0] + "(" + y[1] + ")" for y in sorted(la_paths[k], key=lambda x: x[0])]) + "\t" + s + "\n")
                else:
                    fout.write(
                        k + "\t" + ",".join([y[0] + "(" + y[1] + ")" for y in sorted(ideal_paths[k], key=lambda x: x[0])])
                        + "\t" + "\t" + s + "\n")
                cnt += 1

    print cnt_in, cnt, len(ideal_paths.keys()), len(la_paths.keys())
    print len(only_ideal_lens), len(only_la_lens)
    print "FN ", sorted(only_ideal_lens)[len(only_ideal_lens)/2], sum(only_ideal_lens)*1.0/len(only_ideal_lens), sorted(only_ideal_lens)[len(only_ideal_lens) - 1], len(only_ideal_lens)
    print "FP ", sorted(only_la_lens)[len(only_la_lens)/2], sum(only_la_lens)*1.0/len(only_la_lens), sorted(only_la_lens)[len(only_la_lens) - 1], len(only_la_lens)
    print "Right_num=", right, "Wrong_num=", len(only_ideal_lens) + len(only_la_lens)
    print "FP alone ", fp_alone*1.0/total_alone, fp_alone, total_alone
    precision = right*1.0/(right + len(only_la_lens))
    recall = right*1.0/(right + len(only_ideal_lens))
    precisions.append(precision)
    recalls.append(recall)
    f1scores.append(2*precision * recall/(precision + recall))
    print " threshold", threshold, "precision=", precision, " recall=", recall, " F1-score=", 2* precision * recall/(precision + recall), " tp+fn=", right + len(only_ideal_lens), " tp+fp=", right + len(only_la_lens)
print "Draw Precision.."
plt.figure()
p_line, = plt.plot(thresholds,  precisions, 'r',label = "Precision")
#plt.title('Precision')
#plt.xlabel('edge length threshold')
#plt.ylabel('Another Label')
#plt.savefig(prefix + "_precision.png")

print "Draw Recall.."
#plt.figure()
r_line, = plt.plot(thresholds,  recalls, 'b', label = "Recall")
#plt.title('Recall')
#plt.xlabel('edge length threshold')
#plt.ylabel('Another Label')
#plt.savefig(prefix + "_recall.png")

print "Draw F1score.."
#plt.figure()
f_line, = plt.plot(thresholds,  f1scores, 'g', label = "F-score")
#plt.title('F1-score')
#plt.xlabel('edge length threshold')
#plt.ylabel('Another Label')
plt.legend(handles=[p_line, r_line, f_line])
plt.savefig(prefix + "_multiplescores.png")

# print "Draw FN.."
# print right_0, right, right_small
# plt.figure()
# plt.hist(only_ideal_lens, bins = 50)
# plt.title('FN')
# plt.xlabel('mapping length')
# #plt.ylabel('Another Label')
# plt.savefig(prefix + "_hist_ideal.png")
# print "Draw FP.."
# plt.figure()
# plt.hist(only_la_lens, bins=50)
# plt.title('FP')
# plt.xlabel('mapping length')
# #plt.ylabel('Another Label')
# plt.savefig(prefix + "_hist_la.png")
