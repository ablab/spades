#!/usr/bin/env python
import sys
import json
import unicodedata
import pandas as pd
import datetime
import edlib

import argparse

TRUE_PERFORMANCE = {"ecoli": {"realpb": "88% 87%",
                              "realnp": "86% 87%",
                              "simpb": "99% 83%",
                              "simnp": "96% 88%"}, 
                    "scerevisiae": {"realpb": "56% 87%",
                              "realnp": "61% 83%",
                              "simpb": "98% 83%",
                              "simnp": "85% 83%"},
                    "celegans": {"realpb": "84% 87%",
                              "realnp": "66% 87%",
                              "simpb": "98% 83%",
                              "simnp": "88% 88%"} }

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

    def load_gfa_edges(self, gfa_filename):
        res = {}
        graph = {}
        rev = {"+": "-", "-": "+"}
        with open(gfa_filename, "r") as fin:
            for ln in fin.readlines():
                if ln.startswith("S"):
                    node_id, seq = ln.strip().split("\t")[1:3]
                    res[int(node_id)] = seq
                elif ln.startswith("L"):
                    _, node_id1, pos1, node_id2, pos2, match  = ln.strip().split("\t")
                    if not (node_id1 + pos1 in graph):
                        graph[node_id1 + pos1] = {} 
                    graph[node_id1+pos1][node_id2+pos2] = 1
                    if not (node_id2 + rev[pos2] in graph):
                        graph[node_id2 + rev[pos2]] = {} 
                    graph[node_id2+rev[pos2]][node_id1+rev[pos1]] = 1
        return res, graph

    def load_spaligner_paths(self, filename, reads, stat = "max"):
        res = []
        fin = open(filename, "r")
        for ln in fin.readlines():
            cur_read, seq_starts, seq_ends, e_starts, e_ends, rlen, path, edgelen, seqs = ln.strip().split("\t")
            cur_read = cur_read.split(" ")[0]
            initial_s = [int(x) for x in seq_starts.split(",")] if "," in seq_starts else [int(seq_starts)]
            initial_e = [int(x) for x in seq_ends.split(",")] if "," in seq_ends else [int(seq_ends)]
            mapped_s = [int(x) for x in e_starts.split(",")] if "," in e_starts else [int(e_starts)]
            mapped_e = [int(x) for x in e_ends.split(",")] if "," in e_ends else [int(e_ends)]
            max_ind = []
            for i in xrange(len(initial_s)):
                if stat == "max":
                    if initial_e[i] - initial_s[i] >0.8*len(reads[cur_read]) and (len(max_ind) == 0 or initial_e[i] - initial_s[i] > initial_e[max_ind[0]] - initial_s[max_ind[0]]):
                        max_ind = []
                        max_ind.append(i)
                else:
                    if initial_e[i] - initial_s[i] > min(0.8*len(reads[cur_read]), 1000):
                        max_ind.append(i)
            if len(max_ind) > 0:
                res.append({"r_name": cur_read, "mapping_len": [initial_e[i] - initial_s[i] + 1 for i in max_ind] \
                                   , "s_start": [initial_s[i] for i in max_ind], \
                                     "s_end": [initial_e[i] for i in max_ind], \
                                     "mapped_seq": [seqs.split(";")[i] for i in max_ind]})

        fin.close()
        return res


def save_stats(reads, res_mp):
    df = pd.DataFrame(res_mp)
    df["r"] = df.apply(lambda x: reads[x["r_name"]], axis = 1)
    df["ed"] = df.apply(lambda x: sum([edist([reads[x["r_name"]][x["s_start"][i]: x["s_end"][i]], x["mapped_seq"][i]]) for i in xrange(len(x["mapped_seq"])) ]), axis = 1)
    df["prop_len"] = df.apply(lambda x: sum(x["mapping_len"][i] for i in xrange(len(x["mapping_len"])))*100/len(reads[x["r_name"]]), axis = 1)
    df["prop_ed"] = df.apply(lambda x: x["ed"]*100/sum(x["mapping_len"][i] for i in xrange(len(x["mapping_len"]))), axis = 1)
    res = str(len(df)*100/len(reads))  + "% " + str( 100 - int(df["prop_ed"].mean())) + "%"
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate Alignment Statistics')
    parser.add_argument('-i', '--input', nargs='?', help='Path to folder with all input data', required=True)
    parser.add_argument('-p', '--path', nargs='?', help='Path to folder with results', required=True)
    parser.add_argument('-o', '--orgs', nargs='+', help='Names of organisms to test: ecoli scerevisiae celegans', required=True)

    args = parser.parse_args()
    stat = "max"
    failed = False
    for org in args.orgs:
        res = {}
        for read_type in ["realnp", "realpb", "simpb", "simnp"]:
            reads_file = args.input + "/" + org+ "/input/" + read_type + "2000.fasta"
            dl = DataLoader()
            reads = dl.load_reads(reads_file)
            spaligner_res_file = args.path + "/" + org + "_" + read_type + "/alignment.tsv"
            spaligner_res = dl.load_spaligner_paths(spaligner_res_file, reads, stat)
            res[read_type] = save_stats(reads, spaligner_res)
            if res[read_type] != TRUE_PERFORMANCE[org][read_type]:
                print "Failed: ", org, read_type, res[read_type]
                failed = True
    if failed:
        exit(-1)



