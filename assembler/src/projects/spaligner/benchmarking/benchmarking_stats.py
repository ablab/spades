#!/usr/bin/env python
import sys
import json
import unicodedata
import pandas as pd
import datetime
import edlib

import argparse

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

def make_edit(cur_mp, ideal_seq_offset, graph_seq, ideal_seq):
    if cur_mp["to_length"] == cur_mp["from_length"] and len(cur_mp["sequence"]) == 0:
        graph_seq += ideal_seq[ideal_seq_offset: ideal_seq_offset + int(cur_mp["to_length"])]
    elif cur_mp["to_length"] == cur_mp["from_length"] and len(cur_mp["sequence"]) != 0:
        graph_seq += cur_mp["sequence"]
    elif cur_mp["from_length"] > 0 and cur_mp["to_length"] == 0 and len(cur_mp["sequence"]) != 0:
        graph_seq += cur_mp["sequence"]
    elif cur_mp["from_length"] > 0 and cur_mp["to_length"] == 0:
        print "HAPPENED"
    ideal_seq_offset += int(cur_mp["to_length"])
    return ideal_seq_offset, graph_seq

def make_ga_edit(cur_mp, ideal_seq_offset, graph_seq, ideal_seq):
    if cur_mp["to_length"] == cur_mp["from_length"] and len(cur_mp["sequence"]) == 0:
        graph_seq += ideal_seq[ideal_seq_offset: ideal_seq_offset + int(cur_mp["to_length"])]
    elif cur_mp["to_length"] == cur_mp["from_length"] and len(cur_mp["sequence"]) != 0:
        graph_seq += cur_mp["sequence"]
    elif cur_mp["from_length"] > 0 and cur_mp["to_length"] == 0 and len(cur_mp["sequence"]) != 0:
        graph_seq += cur_mp["sequence"]
    elif cur_mp["from_length"] > 0 and cur_mp["to_length"] == 0:
        print "HAPPENED"
    elif cur_mp["to_length"] != cur_mp["from_length"]:
        graph_seq += cur_mp["sequence"]
    ideal_seq_offset += int(cur_mp["to_length"])
    return ideal_seq_offset, graph_seq

def make_c(s):
    mp = {"A":"T", "T":"A","C":"G", "G":"C"}
    res = "".join([mp[c] for c in list(s)][::-1])
    return res


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

    def load_json_read(self, json_data, ideal_seq, ideal_seq_offset, edges):
        graph_seq = ""
        edge_offset_f = -1
        nodes = []
        for it in json_data["path"]["mapping"]:
            node_id = int(it["position"]["node_id"])
            edge_offset = int(it["position"]["offset"]) if "offset" in it["position"].keys() else 0
            node_seq = edges[node_id] if "is_reverse" not in it["position"].keys() else make_c(edges[node_id])
            edge_offset_f = edge_offset
            ideal_seq_offset_s = ideal_seq_offset
            for mp in it["edit"]:
                cur_mp = {}
                ks = ["from_length", "to_length"]
                for k in ks:
                    cur_mp[k] = mp[k] if k in mp.keys() else "0"
                cur_mp["sequence"] = mp["sequence"] if "sequence" in mp.keys() else ""
                edge_offset_f += int(cur_mp["from_length"])
                ideal_seq_offset, graph_seq = make_ga_edit(cur_mp, ideal_seq_offset, graph_seq, ideal_seq)
            if len(nodes) == 0 or (nodes[-1]["node_id"] != node_id or nodes[-1]["path"][-1]['end'] > edge_offset):
                nodes.append({"node_id": node_id, "node_id_str": str(node_id) + str("+" if "is_reverse" not in it["position"].keys() else "-"),\
                              "path": [{"start": edge_offset, "end": edge_offset_f}], "seq": [{"start": ideal_seq_offset_s, "end": ideal_seq_offset}], "nucs": node_seq})
            else:
                nodes[-1]["path"].append({"start": edge_offset, "end": edge_offset_f})
                nodes[-1]["seq"].append({"start": ideal_seq_offset_s, "end": ideal_seq_offset})
        return nodes

    def extract_max_connected_path(self, nodes, graph):
        prev_node = None
        connected_nodes = []
        cur_nodes = []
        best_len = -1
        cur_len = 0
        for node in nodes:
            if prev_node != None:
                if not prev_node["node_id_str"] in graph or not node["node_id_str"] in graph[prev_node["node_id_str"]]:
                    if cur_len > best_len:
                        best_len = cur_len
                        connected_nodes = cur_nodes
                    cur_len = 0
                    cur_nodes = []
            for p in node["seq"]:
                cur_len += abs(p["end"] - p["start"])
            cur_nodes.append(node)
            prev_node = node

        if cur_len > best_len:
            best_len = cur_len
            connected_nodes = cur_nodes
        return connected_nodes

    def extract_seq(self, connected_nodes):
        graph_seq_full = ""
        for i in xrange(len(connected_nodes)):
            node = connected_nodes[i]
            edge_offset_s = node["path"][0]["start"]
            edge_offset_f = node["path"][-1]["end"]
            graph_seq_full += node["nucs"][edge_offset_s: edge_offset_f]
        return graph_seq_full


    def load_json_paths(self, filename, edges, graph, reads, stat = "max"):
        json_file = open(filename, "r")
        res = []
        res_mp = {}
        for ln in json_file.readlines():
            json_data = json.loads(ln)
            ideal_seq = json_data["sequence"]
            ideal_name = json_data["name"].replace(" ","")
            ideal_seq_offset = 0 if "query_position" not in json_data else int(json_data["query_position"])
            if "path" in json_data.keys():
                nodes = self.load_json_read(json_data, ideal_seq, ideal_seq_offset, edges)
                connected_nodes = self.extract_max_connected_path(nodes, graph)
                graph_seq_full = self.extract_seq(connected_nodes)
                start = connected_nodes[0]["seq"][0]["start"]
                end = connected_nodes[-1]["seq"][-1]["end"]
                if stat == "max":
                    if end - start + 1 > 0.8*len(reads[ideal_name]):
                        if (ideal_name not in res_mp) or (end - start + 1) > res_mp[ideal_name][0]["mapping_len"]:
                            res_mp[ideal_name] = []
                            res_mp[ideal_name].append({"mapping_len": end + 1 - start, "s_start": start, "s_end": end + 1, "mapped_seq": graph_seq_full})
                else:    
                    if end - start + 1 > min(0.8*len(reads[ideal_name]), 1000):
                        if ideal_name not in res_mp:
                            res_mp[ideal_name] = []
                        res_mp[ideal_name].append({"mapping_len": end + 1 - start, "s_start": start, "s_end": end + 1, "mapped_seq": graph_seq_full})
        json_file.close()
        for r in res_mp.keys():
            res.append({"r_name": r, "mapping_len": [x["mapping_len"] for x in res_mp[r]] \
                                   , "s_start": [x["s_start"] for x in res_mp[r]], "s_end": [x["s_end"] for x in res_mp[r]], \
                                     "mapped_seq": [x["mapped_seq"] for x in res_mp[r]]})
        return res

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


def print_stats(reads, res_mp):
    eds = {}
    for name in res_mp.keys():
        print name
        df = pd.DataFrame(res_mp[name]["res"])
        df["r"] = df.apply(lambda x: reads[x["r_name"]], axis = 1)
        df["ed"] = df.apply(lambda x: sum([edist([reads[x["r_name"]][x["s_start"][i]: x["s_end"][i]], x["mapped_seq"][i]]) for i in xrange(len(x["mapped_seq"])) ]), axis = 1)
        df["prop_len"] = df.apply(lambda x: sum(x["mapping_len"][i] for i in xrange(len(x["mapping_len"])))*100/len(reads[x["r_name"]]), axis = 1)
        df["prop_ed"] = df.apply(lambda x: x["ed"]*100/sum(x["mapping_len"][i] for i in xrange(len(x["mapping_len"]))), axis = 1)
        print "Mapped >80%:", len(df)*100/len(reads), "%"
        print "Average Identity:", 100 - int(df["prop_ed"].mean()), "%"
        print "Time:", res_mp[name]["time"]
        print "Memory:", res_mp[name]["memory"]
        print ""
        eds[name] = list(df["prop_ed"])

def get_timedelta(s):
    (h, m, s) = s.split(':')
    d = datetime.timedelta(hours=int(h), minutes=int(m), seconds=int(s))
    return d

def get_time(files_lst):
    res = datetime.timedelta()
    for f in files_lst:
        df = pd.read_table(f)
        df['time'] = df.apply(lambda x: get_timedelta(x['h:m:s']), axis=1)
        res += min(df['time'])
    return res

def get_memory(files_lst):
    res = 0
    for f in files_lst:
        df = pd.read_table(f)
        res = max(res,max(df['max_rss']))
    return res

def save_fasta(aligner_res, filename):
    with open(filename, "w") as fout:
        for i in xrange(len(aligner_res)):
            fout.write(">" + aligner_res[i]["r_name"] + "\n" + aligner_res[i]["mapped_seq"] + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate Alignment Statistics')
    parser.add_argument('-p', '--path', nargs='?', help='Path to folder with results for all aligners', required=True)
    parser.add_argument('-a', '--aligners', nargs='+', help='Names of aligners to test: SPAligner vg_xdrop vg_ordinary GraphAligner', required=True)
    parser.add_argument('-o', '--orgs', nargs='+', help='Names of organisms to test: ecoli scerevisiae celegans', required=True)

    args = parser.parse_args()
    datapath = args.path
    aligners = {"SPAligner": 0, "vg_xdrop": 0, "vg_ordinary": 0, "GraphAligner":0}
    for al in args.aligners:
        aligners[al] = 1
    stat = "max"
    for org in args.orgs:
        for read_type in ["realnp2000", "realpb2000", "simpb2000", "simnp2000"]:
            org_path = datapath + "/" + org + "/"
            reads_file = org_path + "input/" + read_type + ".fasta"
            print org, read_type
            dl = DataLoader()
            reads = dl.load_reads(reads_file)
            mp = {}
            for al in aligners.keys():
                if aligners[al] == 1:
                    if al.startswith("SPAligner"):
                        spaligner_res_file = org_path + al + "/output/aln_" + read_type + ".tsv"
                        spaligner_res = dl.load_spaligner_paths(spaligner_res_file, reads, stat)
                        log_files = [org_path + al + "/benchmark/align_" + read_type + ".tsv"]
                        mp[al] = {"time": get_time(log_files), "memory": get_memory(log_files)}
                        mp[al]["res"] = spaligner_res

                    if al.startswith("GraphAligner"):
                        graphaligner_edges_gfa = org_path + al + "/tmp/graph_idfix.gfa"
                        graphaligner_res_file = org_path + al + "/output/aln_" + read_type + "_selected.json"
                        [graphaligner_edges, graphaligner_graph] = dl.load_gfa_edges(graphaligner_edges_gfa)
                        graphaligner_res = dl.load_json_paths(graphaligner_res_file, graphaligner_edges, graphaligner_graph, reads, stat)
                        log_files = [org_path + al + "/benchmark/mummer_pipe_" + read_type + ".tsv",\
                                    org_path + al + "/benchmark/align_" + read_type + ".tsv",\
                                    org_path + al + "/benchmark/posprocess_" + read_type + ".tsv"]
                        mp[al] = {"time": get_time(log_files), "memory": get_memory(log_files)}
                        mp[al]["res"] = graphaligner_res

                    if al.startswith("vg_xdrop"):
                        vg_edges_gfa = org_path + al + "/tmp/graph.split.gfa"
                        vg_res_file = org_path + al + "/output/aln_" + read_type + "_xdrop.json"
                        [vg_edges, vg_graph] = dl.load_gfa_edges(vg_edges_gfa)
                        vg_res = dl.load_json_paths(vg_res_file, vg_edges, vg_graph, reads, stat)
                        log_files = [org_path + al + "/benchmark/build_index.tsv",\
                                    org_path + al + "/benchmark/mapping_" + read_type + "_xdrop.tsv"]
                        if org == "celegans":
                            log_files.append(org_path + al + "/benchmark/prune_graph.tsv")
                        mp[al] = {"time": get_time(log_files), "memory": get_memory(log_files)}
                        mp[al]["res"]  = vg_res

                    if al.startswith("vg_ordinary"):
                        vg_edges_gfa = org_path + al + "/tmp/graph.split.gfa"
                        vg_res_file = org_path + al + "/output/aln_" + read_type + "_ordinary.json"
                        [vg_edges, vg_graph] = dl.load_gfa_edges(vg_edges_gfa)
                        vg_res = dl.load_json_paths(vg_res_file, vg_edges, vg_graph, reads, stat)
                        log_files = [org_path + al + "/benchmark/build_index.tsv",\
                                    org_path + al + "/benchmark/mapping_" + read_type + "_ordinary.tsv"]
                        if org == "celegans":
                            log_files.append(org_path + al + "/benchmark/prune_graph.tsv")
                        mp[al] = {"time": get_time(log_files), "memory": get_memory(log_files)}
                        mp[al]["res"]  = vg_res
            print_stats(reads, mp)


