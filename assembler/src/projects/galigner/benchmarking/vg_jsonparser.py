import json
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

def load_edges(gfa_filename):
    res = {}
    with open(gfa_filename, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith("S"):
                _, node_id, seq = ln.strip().split("\t")
                res[int(node_id)] = seq
    return res

def load_reads(filename):
    res = {}
    key = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip()
                res[key] = ""
            else:
                res[key] += ln.strip()
    return res

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


filename = "/home/tdvorkina/soft/vg/ecoli_real/filtered_subreads_mapped_10k.fasta.json"
gfa_filename = "/home/tdvorkina/soft/vg/ecoli_real/assembly_graph_with_scaffolds.split.gfa"

reads = load_reads("/home/tdvorkina/soft/vg/ecoli_real/filtered_subreads_mapped_10k.fasta")
edges = load_edges(gfa_filename)

json_file = open(filename, "r")
vg_ed = {}
vg_ed_lst = []
vg_ed_per_lst = []
for ln in json_file.readlines():
    json_data = json.loads(ln)
    ideal_seq = json_data["sequence"]
    ideal_name = json_data["name"]
    ideal_seq_offset = 0
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
            ideal_seq_offset, graph_seq = make_edit(cur_mp, ideal_seq_offset, graph_seq, ideal_seq)
        if len(nodes) == 0 or nodes[-1]["node_id"] != node_id:
            nodes.append({"node_id": node_id, "path": [{"start": edge_offset, "end": edge_offset_f}], "seq": [{"start": ideal_seq_offset_s, "end": ideal_seq_offset}], "nucs": node_seq})
        else:
            nodes[-1]["path"].append({"start": edge_offset, "end": edge_offset_f})
            nodes[-1]["seq"].append({"start": ideal_seq_offset_s, "end": ideal_seq_offset})

    graph_seq_full = ""
#    print len(ideal_seq), ideal_seq_offset
    for i in xrange(len(nodes)):
        node = nodes[i]
#        print "edge_id=", node["node_id"], "  path=", "; ".join(str(x["start"]) + "-" + str(x["end"]) for x in node["path"])
        edge_offset_s = 0
        edge_offset_f = len(node["nucs"])
        if i == 0:
            edge_offset_s = node["path"][0]["start"]
        if i == len(nodes) - 1:
            edge_offset_f = node["path"][-1]["end"]
        graph_seq_full += node["nucs"][edge_offset_s: edge_offset_f + 1]

    vg_ed[ideal_name] = edist([ideal_seq, graph_seq_full])*100/len(ideal_seq)
    vg_ed_lst.append(edist([ideal_seq, graph_seq_full]))
    vg_ed_per_lst.append(edist([ideal_seq, graph_seq_full])*100/len(ideal_seq))
    # if edist([ideal_seq, graph_seq])*100/len(ideal_seq) < edist([ideal_seq, graph_seq_full])*100/len(ideal_seq):
    #     print edist([ideal_seq, graph_seq])*100/len(ideal_seq), edist([ideal_seq, graph_seq_full])*100/len(ideal_seq),  " identity=", json_data["identity"], " score=", json_data["score"] 
#    break

galigner_tsv = "/home/tdvorkina/results//vg_comparison/master_2018-05-18_11-48-24_E.coli_real_synth_dijkstra_bwa200.tsv"

galigner_ed = {}
galigner_ed_lst = []
galigner_ed_per_lst = []
with open(galigner_tsv, "r") as fin:
    for ln in fin.readlines():
        cur_read, seq_start, seq_end, rlen, path_dirty, edgelen, ss, ed = ln.strip().split("\t")
        ed = edist([reads[cur_read], ss])
        galigner_ed[cur_read] = int(ed)*100/int(rlen)
        galigner_ed_lst.append(int(ed))
        galigner_ed_per_lst.append(int(ed)*100/int(rlen))

print "Mapped num galigner=", len(galigner_ed), " vg=", len(vg_ed)
print "Median ed galigner=", sorted([galigner_ed[x] for x in galigner_ed.keys()])[len(galigner_ed)/2], " vg=", sorted([vg_ed[x] for x in vg_ed.keys()])[len(vg_ed)/2]

plt.hist(vg_ed_lst, alpha=0.5, label="vg", bins=100)
plt.hist(galigner_ed_lst, alpha=0.5, label="galigner", bins=100)
plt.legend(loc='upper right')
plt.title("Edit distance")
plt.savefig("vg_ga_ed_real_ecoli.png")

plt.figure()
plt.hist(vg_ed_per_lst, alpha=0.5, label="vg", bins=100)
plt.hist(galigner_ed_per_lst, alpha=0.5, label="galigner", bins=100)
plt.legend(loc='upper right')
plt.title("Edit distance / read length (%)")
plt.savefig("vg_ga_ed_per_real_ecoli.png")