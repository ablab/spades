import sys
import json
import unicodedata
import pandas as pd
import edlib
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def edist(lst):
    #return editdistance.eval(lst[0], lst[1])
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

def edist_locations(lst):
    #return editdistance.eval(lst[0], lst[1])
    result = edlib.align(str(lst[0]), str(lst[1]), mode="HW", task="locations", additionalEqualities=[('U', 'T')
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
    return result["editDistance"], result["locations"][0][0], result["locations"][0][1] 

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

    def load_vg_edges(self, gfa_filename):
        res = {}
        graph = {}
        with open(gfa_filename, "r") as fin:
            for ln in fin.readlines():
                if ln.startswith("S"):
                    _, node_id, seq = ln.strip().split("\t")
                    res[int(node_id)] = seq
                elif ln.startswith("L"):
                    _, node_id1, pos1, node_id2, pos2, match  = ln.strip().split("\t")
                    if not (node_id1 + pos1 in graph):
                        graph[node_id1 + pos1] = {} 
                    graph[node_id1+pos1][node_id2+pos2] = 1

        return res, graph

    def load_vg_paths(self, filename, edges, graph, reads):
        json_file = open(filename, "r")
        res = []
        res_mp = {}
        for ln in json_file.readlines():
            json_data = json.loads(ln)
            ideal_seq = json_data["sequence"]
            ideal_name = json_data["name"].replace(" ","")
            ideal_seq_offset = 0
            graph_seq = ""
            edge_offset_f = -1
            nodes = []
            if "path" in json_data.keys():
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
                        ideal_seq_offset, graph_seq = make_edit(cur_mp, ideal_seq_offset, graph_seq, ideal_seq)
                    if len(nodes) == 0 or nodes[-1]["node_id"] != node_id:
                        nodes.append({"node_id": node_id, "node_id_str": str(node_id) + str("+" if "is_reverse" not in it["position"].keys() else "-"),\
                                      "path": [{"start": edge_offset, "end": edge_offset_f}], "seq": [{"start": ideal_seq_offset_s, "end": ideal_seq_offset}], "nucs": node_seq})
                    else:
                        nodes[-1]["path"].append({"start": edge_offset, "end": edge_offset_f})
                        nodes[-1]["seq"].append({"start": ideal_seq_offset_s, "end": ideal_seq_offset})

                prev_node = None
                connected_nodes = []
                cur_nodes = []
                best_len = -1
                cur_len = 0
                for node in nodes:
                    if prev_node != None:
                        if not prev_node["node_id_str"] in graph or not node["node_id_str"] in graph[prev_node["node_id_str"]]:
                            if ideal_name == "S1_26776":
                                print "Disconnected: ", ideal_name, cur_len, best_len
                                print cur_nodes 
                            if cur_len > best_len:
                                best_len = cur_len
                                connected_nodes = cur_nodes
                            cur_len = 0
                            cur_nodes = []
                    for p in node["seq"]:
                        cur_len += abs(p["end"] - p["start"])
                    cur_nodes.append(node)
                    #print cur_len
                    prev_node = node

                if cur_len > best_len:
                    best_len = cur_len
                    connected_nodes = cur_nodes
                graph_seq_full = ""
                graph_seq = ""
                #print ideal_name, best_len, len(ideal_seq), len(reads[ideal_name])
                #print len(connected_nodes)
                for i in xrange(len(connected_nodes)):
                    node = connected_nodes[i]
                    edge_offset_s = 0
                    edge_offset_f = len(node["nucs"])
                    if i == 0:
                        edge_offset_s = node["path"][0]["start"]
                    if i == len(nodes) - 1:
                        edge_offset_f = node["path"][-1]["end"]
                    #graph_seq_full = edges[node["node_id"]] if node["node_id_str"].endswith("+") else make_c(edges[node["node_id"]])
                    graph_seq_full += node["nucs"][edge_offset_s: edge_offset_f + 1]
                    #print edge_offset_s, edge_offset_f + 1, len(node["nucs"]), node["node_id_str"], len(edges[node["node_id"]])
                    # for p in node["seq"]:
                    #     graph_seq += node["nucs"][p["start"]:p["end"]]

                start = connected_nodes[0]["seq"][0]["start"]
                end = connected_nodes[-1]["seq"][-1]["end"]
    
                if end - start + 1 > 0.*len(reads[ideal_name]):
                    if (ideal_name not in res_mp) or (end - start + 1) > res_mp[ideal_name]["mapping_len"]:
                        #print ideal_name, end - start + 1, len(reads[ideal_name])
                        res_mp[ideal_name] = {"mapping_len": end + 1 - start, "s_start": start, "s_end": end + 1, "mapped_seq": graph_seq_full}
        json_file.close()
        for r in res_mp.keys():
            res.append({"r_name": r, "mapping_len": res_mp[r]["mapping_len"], "s_start": res_mp[r]["s_start"], "s_end": res_mp[r]["s_end"], \
                                     "mapped_seq": res_mp[r]["mapped_seq"]})
        return res

    def load_galigner_paths(self, filename, reads):
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
            if initial_e[max_ind] - initial_s[max_ind] + 1 > 0.8*len(reads[cur_read]):
                res.append({"r_name": cur_read, "mapping_len": initial_e[max_ind] - initial_s[max_ind] + 1, \
                                "s_start": initial_s[max_ind], "s_end": initial_e[max_ind], \
                                "mapped_seq": seqs.split(";")[max_ind]})
            else:
                print cur_read
                print initial_e[max_ind] - initial_s[max_ind], len(reads[cur_read])
        fin.close()
        return res


def print_stats(reads, res_mp):
    eds = {}
    for name in res_mp.keys():
        print name
        df = pd.DataFrame(res_mp[name])
        df["r"] = df.apply(lambda x: reads[x["r_name"]], axis = 1)
        df["ed"] = df.apply(lambda x: edist([reads[x["r_name"]][x["s_start"]: x["s_end"]], x["mapped_seq"]]), axis = 1)
        df["prop_len"] = df.apply(lambda x: x["mapping_len"]*100/len(reads[x["r_name"]]), axis = 1)
        df["prop_ed"] = df.apply(lambda x: x["ed"]*100/x["mapping_len"], axis = 1)
        print "Mapped:", len(df)*100/len(reads), "%"
        print "Mean len:", int(df["prop_len"].mean()), "%"
        print "Mean ed:", int(df["prop_ed"].mean()), "%"
        print "Median len:", int(df["prop_len"].median()), "%"
        print "Median ed:", int(df["prop_ed"].median()), "%"
        print ""
        eds[name] = list(df["prop_ed"])

    for k in eds.keys():
        plt.hist(eds[k], alpha=0.5, label=k, bins=20)
    plt.legend(loc='upper right')
    plt.title("Edit distance / read length (%)")
    plt.savefig("ed_per.png")

def save_fasta(aligner_res, filename):
    with open(filename, "w") as fout:
        for i in xrange(len(aligner_res)):
            fout.write(">" + aligner_res[i]["r_name"] + "\n" + aligner_res[i]["mapped_seq"] + "\n")


if __name__ == "__main__":
    # reads_file = "/Sid/tdvorkina/gralign/C.elegans_synth/benchmarking/sim_pacbio/sim_pacbio_len500_100.fasta"
    # galigner_res_file = "/home/tdvorkina/results//benchmarking_test/run2_2018-07-03_17-43-17_C.elegans_simpb100_dijkstra_bwa0_ends_inf_extended.tsv"
    # edges_gfa = "/home/tdvorkina/soft/vg/celegans_sim_pacbio/assembly_graph_with_scaffolds_wp.split.gfa"
    # vg_res_file = "/home/tdvorkina/soft/vg/celegans_sim_pacbio/sim_pacbio_len500_100.json"

    reads_file = "/Sid/tdvorkina/gralign/E.coli_synth/benchmarking/sim_pacbio/sim_pacbio_len500_100.fasta"
    galigner_res_file = "/home/tdvorkina/results//benchmarking_test/run2_2018-07-03_17-43-17_E.coli_simpb100_dijkstra_bwa0_ends_inf.tsv"
    #edges_gfa = "/home/tdvorkina/soft/vg/ecoli_sim_pacbio/assembly_graph_with_scaffolds_wp.split.gfa"
    #vg_res_file = "/home/tdvorkina/soft/vg/ecoli_sim_pacbio/sim_pacbio_len500_100.json"
    graphaligner_edges_gfa = "/home/tdvorkina/soft/GraphAligner/tmp/graph_idfix.gfa"
    graphaligner_res_file = "/home/tdvorkina/soft/GraphAligner/output/aln_sim_pacbio_len500_100_graph_all.json"

    # reads_file = "/Sid/tdvorkina/gralign/E.coli_synth/benchmarking/real_pacbio/real_pacbio_len500_100.fasta"
    # galigner_res_file =  "/home/tdvorkina/results//benchmarking_test/run2_2018-07-03_17-43-17_E.coli_realpb100_dijkstra_bwa0_ends_inf.tsv"
    # edges_gfa = "/home/tdvorkina/soft/vg/ecoli_sim_pacbio/assembly_graph_with_scaffolds_wp.split.gfa"
    # vg_res_file = "/home/tdvorkina/soft/vg/ecoli_real_pacbio/real_pacbio_len500_100.json"

    # reads_file = "/Sid/tdvorkina/gralign/C.elegans_synth/benchmarking/real_pacbio/real_pacbio_len500_100.fasta"
    # galigner_res_file = "/home/tdvorkina/results//benchmarking_test/run2_2018-07-03_17-43-17_C.elegans_realpb100_dijkstra_bwa0_ends_inf.tsv"
    # edges_gfa = "/home/tdvorkina/soft/vg/celegans_sim_pacbio/assembly_graph_with_scaffolds_wp.split.gfa"
    # vg_res_file = "/home/tdvorkina/soft/vg/celegans_real_pacbio/real_pacbio_len500_100.json"
    
    dl = DataLoader()
    reads = dl.load_reads(reads_file)
    galigner_res = dl.load_galigner_paths(galigner_res_file, reads)
    # [vg_edges, vg_graph] = dl.load_vg_edges(edges_gfa)
    # vg_res = dl.load_vg_paths(vg_res_file, vg_edges, vg_graph, reads)
    [graphaligner_edges, graphaligner_graph] = dl.load_vg_edges(graphaligner_edges_gfa)
    graphaligner_res = dl.load_vg_paths(graphaligner_res_file, graphaligner_edges, graphaligner_graph, reads)
    print_stats(reads, {"GAligner": galigner_res, "GraphAligner": graphaligner_res})
    filename = "/".join(vg_res_file.split("/")[:-1]) + "/alignments.fasta"
    print "Saving to", filename
    save_fasta(vg_res, filename)
    print "Draw reads len distribution.."
    lens = [len(reads[k]) for k in reads.keys()]
    plt.hist(lens)
    plt.xlabel("Read length")
    plt.savefig("./len_dist.png")