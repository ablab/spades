#!/usr/bin/env python2
import sys

def load_lineage(file):
    fin = open(file, "r")
    res = {"superkingdom":set(), "phylum":set(), "class":set(), "order":set(), "family":set(), "genus": set(), "species": set()}
    num = 0
    for ln in fin.readlines():
        num += 1
        lst = ln.strip().split("\t")[1:]
        for it in lst:
            #print it
            tp, name = it.split(":")
            if tp not in res.keys():
                print "WARNING ", tp 
                continue
            res[tp].add(name)
    fin.close()
    return res

def load_lineage_lst(file):
    fin = open(file, "r")
    res = []
    for ln in fin.readlines():
        if len(ln.strip().split("\t")[1:]) > 0:
            res.append(ln.strip().split("\t")[1:])
    fin.close()
    return res

file1 = sys.argv[1]
file2 = sys.argv[2]

lineage1 = load_lineage(file1)
lineage2 = load_lineage(file2)

tps = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
for tp in tps:
    print tp, " ideal num ",  len(lineage1[tp]), " mapped num ",len(lineage2[tp]), " intersect num",len(lineage1[tp].intersection(lineage2[tp]))
    #print sorted(lineage1[tp] - lineage2[tp])
    #print sorted(lineage1[tp].intersection(lineage2[tp]))

l1 = load_lineage_lst(file1)
l2 = load_lineage_lst(file2)

for i in xrange(len(tps)):
    tp = tps[i]
    true = set()
    maybe = 0
    s1 = set()
    s2 = set()
    for it in l1:
        ind = 0
        #while not it[ind].startswith(""):
        str1 = " ".join(it[:i + 1])
        s1.add(str1)
        #print str1
        for it2 in l2:
            if len(it2) < i + 1:
               continue
            str2 = " ".join(it2[:i + 1])    
            s2.add(str2)
            if str1 == str2:
               true.add(str1)

    for it2 in l2:
        if len(it2) < i + 1:
            maybe += 1
    print tps[i], "ideal num ", len(s1), " Mapped num ", len(s2), " Intersec ", len(true), " Maybe ", maybe
